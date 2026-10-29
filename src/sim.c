#include "sim.h"
#include "util.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

/* ───────────────────────────────────────────────────────── lifecycle ── */

int sim_init(Sim *s, const Params *p)
{
    memset(s, 0, sizeof *s);
    s->p = *p;
    s->cols = p->cols; s->rows = p->rows; s->n = p->cols * p->rows;
    s->held = -1;
    s->rng = (uint64_t)p->seed * 0x2545F4914F6CDD1Dull + 1;

    size_t n = (size_t)s->n;
    double **arrs[] = { &s->ux, &s->uy, &s->vx, &s->vy, &s->ax, &s->ay,
                        &s->m, &s->inv_m, &s->gamma, &s->decay, &s->anchor_k };
    for (size_t i = 0; i < sizeof arrs / sizeof *arrs; i++)
        if (!(*arrs[i] = calloc(n, sizeof(double)))) return -1;
    if (!(s->fixed = calloc(n, 1))) return -1;
    s->spr_cap = s->n * 4 + 8;
    if (!(s->spr = calloc((size_t)s->spr_cap, sizeof *s->spr))) return -1;

    sim_set_masses(s);
    sim_rebuild(s);
    return 0;
}

void sim_free(Sim *s)
{
    free(s->ux); free(s->uy); free(s->vx); free(s->vy); free(s->ax); free(s->ay);
    free(s->m); free(s->inv_m); free(s->gamma); free(s->decay); free(s->anchor_k);
    free(s->fixed); free(s->spr);
    memset(s, 0, sizeof *s);
}

int sim_defect_index(const Sim *s)
{
    if (s->p.disorder != DIS_LIGHT && s->p.disorder != DIS_HEAVY) return -1;
    return sim_idx(s, s->cols / 2, s->rows / 2);
}

void sim_set_masses(Sim *s)
{
    const Params *p = &s->p;
    for (int i = 0; i < s->n; i++) {
        double m = p->mass;
        if (p->disorder == DIS_RANDOM) m *= 1.0 + p->disorder_width * rng_sym(&s->rng);
        s->m[i] = m;
    }
    int d = sim_defect_index(s);
    if (d >= 0) s->m[d] = p->mass * (p->disorder == DIS_LIGHT ? p->defect_light : p->defect_heavy);
    for (int i = 0; i < s->n; i++) s->inv_m[i] = 1.0 / s->m[i];
    sim_update_damping(s);
}

double sim_mean_mass(const Sim *s)
{
    double sum = 0; for (int i = 0; i < s->n; i++) sum += s->m[i];
    return sum / s->n;
}

double sim_min_mass(const Sim *s)
{
    double mn = s->m[0]; for (int i = 1; i < s->n; i++) if (s->m[i] < mn) mn = s->m[i];
    return mn;
}

double sim_omega0(const Sim *s) { return sqrt(s->p.k / s->p.mass); }

/* ───────────────────────────────────────────────────────── topology ── */

static void add_spring(Sim *s, int a, int b, double rx, double ry, double k, int wrap)
{
    double len = sqrt(rx * rx + ry * ry);
    s->spr[s->nspr++] = (Spring){ a, b, rx, ry, k, len * (1.0 - s->prestrain_eff), wrap };
}

void sim_rebuild(Sim *s)
{
    const Params *p = &s->p;
    const int C = s->cols, R = s->rows;
    const int per = p->bc == BC_PERIODIC;

    s->prestrain_eff = (p->bc == BC_CLAMPED || p->bc == BC_PERIODIC) ? p->prestrain : 0.0;
    s->nspr = 0;

    for (int r = 0; r < R; r++)
        for (int c = 0; c < C; c++) {
            int a = sim_idx(s, c, r);
            int cr = c + 1, cl = c - 1, rd = r + 1;
            int wr = cr >= C, wl = cl < 0, wd = rd >= R;
            if (per) { cr %= C; cl = (cl + C) % C; rd %= R; }
            int okr = per || !wr, okl = per || !wl, okd = per || !wd;

            if (okr)       add_spring(s, a, sim_idx(s, cr, r),  1, 0, p->k, wr);
            if (okd)       add_spring(s, a, sim_idx(s, c, rd),  0, 1, p->k, wd);
            if (p->diagonal) {
                if (okr && okd) add_spring(s, a, sim_idx(s, cr, rd),  1, 1, p->kdiag, wr || wd);
                if (okl && okd) add_spring(s, a, sim_idx(s, cl, rd), -1, 1, p->kdiag, wl || wd);
            }
        }

    /* boundary nodes */
    for (int r = 0; r < R; r++)
        for (int c = 0; c < C; c++) {
            int i = sim_idx(s, c, r);
            int edges = (c == 0) + (c == C - 1) + (r == 0) + (r == R - 1);
            s->fixed[i]    = (p->bc == BC_CLAMPED && edges > 0);
            s->anchor_k[i] = (p->bc == BC_PINNED) ? p->k * edges : 0.0;
            if (s->fixed[i]) s->ux[i] = s->uy[i] = s->vx[i] = s->vy[i] = 0.0;
        }
    if (s->held >= 0 && s->fixed[s->held]) s->held = -1;

    /* reference-state potential energy (nonzero under prestress) */
    double pe = 0;
    for (int k = 0; k < s->nspr; k++) {
        double len = sqrt(s->spr[k].rx * s->spr[k].rx + s->spr[k].ry * s->spr[k].ry);
        double d = len - s->spr[k].L0;
        pe += 0.5 * s->spr[k].k * d * d;
    }
    s->pe_ref = pe;

    /* Gershgorin bound on omega_max^2 */
    double *row = calloc((size_t)s->n, sizeof *row);
    if (row) {
        for (int k = 0; k < s->nspr; k++) {
            const Spring *sp = &s->spr[k];
            double len = sqrt(sp->rx * sp->rx + sp->ry * sp->ry);
            double t = sp->k * (len - sp->L0) / len;
            double w = 2.0 * (sp->k + fabs(t));
            row[sp->a] += w; row[sp->b] += w;
        }
        double mx = 0;
        for (int i = 0; i < s->n; i++) {
            double v = (row[i] + s->anchor_k[i]) * s->inv_m[i];
            if (v > mx) mx = v;
        }
        s->omega_max = sqrt(mx);
        free(row);
    }

    sim_update_damping(s);
    sim_compute_forces(s);
}

void sim_update_damping(Sim *s)
{
    const Params *p = &s->p;
    const int C = s->cols, R = s->rows, W = p->absorb_width;
    for (int r = 0; r < R; r++)
        for (int c = 0; c < C; c++) {
            int i = sim_idx(s, c, r);
            double g = p->damping;
            if (p->absorb) {
                int d = c;
                if (r < d) d = r;
                if (C - 1 - c < d) d = C - 1 - c;
                if (R - 1 - r < d) d = R - 1 - r;
                if (d < W) { double x = (double)(W - d) / W; g += p->absorb_gamma * x * x; }
            }
            s->gamma[i] = g;
            s->decay[i] = exp(-g * p->dt * 0.5 * s->inv_m[i]);
        }
}

void sim_reset(Sim *s)
{
    memset(s->ux, 0, (size_t)s->n * sizeof(double));
    memset(s->uy, 0, (size_t)s->n * sizeof(double));
    memset(s->vx, 0, (size_t)s->n * sizeof(double));
    memset(s->vy, 0, (size_t)s->n * sizeof(double));
    s->held = -1;
    s->t = 0; s->steps = 0;
    sim_compute_forces(s);
}

/* ───────────────────────────────────────────────────────── dynamics ── */

void sim_compute_forces(Sim *s)
{
    const int n = s->n;
    double *ux = s->ux, *uy = s->uy, *ax = s->ax, *ay = s->ay;

    for (int i = 0; i < n; i++) {
        ax[i] = -s->anchor_k[i] * ux[i];
        ay[i] = -s->anchor_k[i] * uy[i];
    }
    for (int k = 0; k < s->nspr; k++) {
        const Spring *sp = &s->spr[k];
        int a = sp->a, b = sp->b;
        double dx = sp->rx + ux[b] - ux[a];
        double dy = sp->ry + uy[b] - uy[a];
        double dist = sqrt(dx * dx + dy * dy);
        if (dist < 1e-12) continue;
        double f = sp->k * (dist - sp->L0) / dist;
        double fx = f * dx, fy = f * dy;
        ax[a] += fx; ay[a] += fy;
        ax[b] -= fx; ay[b] -= fy;
    }
    for (int i = 0; i < n; i++) { ax[i] *= s->inv_m[i]; ay[i] *= s->inv_m[i]; }
}

static void damp_half(Sim *s)
{
    for (int i = 0; i < s->n; i++) { s->vx[i] *= s->decay[i]; s->vy[i] *= s->decay[i]; }
}

void sim_step(Sim *s)
{
    const double dt = s->p.dt, h = 0.5 * dt;
    const int n = s->n, held = s->held;

    damp_half(s);
    for (int i = 0; i < n; i++) {
        if (s->fixed[i] || i == held) continue;
        s->vx[i] += h * s->ax[i];  s->vy[i] += h * s->ay[i];
        s->ux[i] += dt * s->vx[i]; s->uy[i] += dt * s->vy[i];
    }
    if (held >= 0) {
        s->ux[held] = s->held_ux; s->uy[held] = s->held_uy;
        s->vx[held] = s->held_vx; s->vy[held] = s->held_vy;
    }
    sim_compute_forces(s);
    for (int i = 0; i < n; i++) {
        if (s->fixed[i] || i == held) continue;
        s->vx[i] += h * s->ax[i]; s->vy[i] += h * s->ay[i];
    }
    damp_half(s);
    if (held >= 0) { s->vx[held] = s->held_vx; s->vy[held] = s->held_vy; }

    s->steps++;
    s->t = s->steps * dt;
}

void sim_energy(const Sim *s, double *ke, double *pe)
{
    double k = 0, v = 0;
    for (int i = 0; i < s->n; i++) {
        k += 0.5 * s->m[i] * (s->vx[i] * s->vx[i] + s->vy[i] * s->vy[i]);
        v += 0.5 * s->anchor_k[i] * (s->ux[i] * s->ux[i] + s->uy[i] * s->uy[i]);
    }
    for (int j = 0; j < s->nspr; j++) {
        const Spring *sp = &s->spr[j];
        double dx = sp->rx + s->ux[sp->b] - s->ux[sp->a];
        double dy = sp->ry + s->uy[sp->b] - s->uy[sp->a];
        double d = sqrt(dx * dx + dy * dy) - sp->L0;
        v += 0.5 * sp->k * d * d;
    }
    *ke = k; *pe = v - s->pe_ref;
}

void sim_lattice(const Sim *s, Lattice *L)
{
    double e = s->prestrain_eff;
    L->k = s->p.k; L->kd = s->p.kdiag; L->diag = s->p.diagonal;
    L->tc = s->p.k * e;          /* T/|r| = k(|r| - L0)/|r| = k e */
    L->td = s->p.kdiag * e;
    L->m = sim_mean_mass(s);
}

/* ───────────────────────────────────────────────────────── perturbations ── */

void sim_kick_random(Sim *s, double ampx, double ampy)
{
    for (int i = 0; i < s->n; i++) {
        if (s->fixed[i] || i == s->held) continue;
        s->vx[i] += ampx * rng_sym(&s->rng);
        s->vy[i] += ampy * rng_sym(&s->rng);
    }
}

void sim_displace(Sim *s, int i, double dx, double dy)
{
    if (i < 0 || i >= s->n || s->fixed[i]) return;
    s->ux[i] += dx; s->uy[i] += dy;
    s->vx[i] = s->vy[i] = 0.0;
    sim_compute_forces(s);
}

int sim_nearest(const Sim *s, double x, double y, double maxd)
{
    int best = -1; double best2 = maxd * maxd;
    for (int r = 0; r < s->rows; r++)
        for (int c = 0; c < s->cols; c++) {
            int i = sim_idx(s, c, r);
            if (s->fixed[i]) continue;
            double dx = c + s->ux[i] - x, dy = r + s->uy[i] - y, d2 = dx * dx + dy * dy;
            if (d2 < best2) { best2 = d2; best = i; }
        }
    return best;
}

void sim_excite_mode(Sim *s, int branch, double strain)
{
    const Boundary bc = s->p.bc;
    const int C = s->cols, R = s->rows;
    int wx = line_w0(bc) == 0 ? 1 : line_w0(bc);     /* lowest nonzero x mode */
    int wy = line_w0(bc);                            /* lowest y profile      */
    double qx = line_q(bc, C, wx), qy = line_q(bc, R, wy);

    Lattice L; sim_lattice(s, &L);
    double w[2], e[2][2];
    branches(&L, qx, qy, w, e);

    double fx_max = 0, fy_max = 0;
    for (int c = 0; c < C; c++) { double v = fabs(line_phi(bc, C, wx, 0, c)); if (v > fx_max) fx_max = v; }
    for (int r = 0; r < R; r++) { double v = fabs(line_phi(bc, R, wy, 0, r)); if (v > fy_max) fy_max = v; }

    double A = strain / qx;
    for (int r = 0; r < R; r++)
        for (int c = 0; c < C; c++) {
            int i = sim_idx(s, c, r);
            double f = A * line_phi(bc, C, wx, 0, c) / fx_max * line_phi(bc, R, wy, 0, r) / fy_max;
            s->ux[i] = f * e[branch][0];
            s->uy[i] = f * e[branch][1];
            s->vx[i] = s->vy[i] = 0.0;
            if (s->fixed[i]) s->ux[i] = s->uy[i] = 0.0;
        }
    s->held = -1;
    sim_compute_forces(s);
}
