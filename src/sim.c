#include "sim.h"
#include "util.h"

#include <math.h>
#include <stdio.h>
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
                        &s->m, &s->inv_m, &s->gamma, &s->decay, &s->anchor_k,
                        &s->oux, &s->ouy, &s->sx, &s->sy, &s->svx, &s->svy };
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
    free(s->oux); free(s->ouy); free(s->sx); free(s->sy); free(s->svx); free(s->svy);
    free(s->fixed); free(s->spr); free(s->chain); free(s->chain_off); free(s->tri); free(s->vfac);
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

#define ROD_TOL      1e-8      /* SHAKE: relative |d^2 - L^2| / L^2      */
#define ROD_MAXIT    30         /* per step; undisturbed steps need ~3-5 */
#define PROJ_MAXIT   2000       /* one-off projections */
#define PROJ_TOL     1e-7      /* a projection that ends above this failed */
#define BLOCK_TOL    1e-4      /* a step that ends above this is undone    */


static void add_spring(Sim *s, int a, int b, double rx, double ry, double k, int wrap)
{
    double len = sqrt(rx * rx + ry * ry);
    int diag = rx != 0 && ry != 0;
    int rod = diag ? s->p.rod_diag : s->p.rod_card;
    /* rods keep the reference length; prestrain applies to springs only */
    double L0 = rod ? len : len * (1.0 - s->prestrain_eff);
    s->spr[s->nspr++] = (Spring){ a, b, rx, ry, k, L0, wrap, rod, 0.0 };
    if (rod) s->nrods++;
}

static inline double wgt(const Sim *s, int i)
{
    return (s->fixed[i] || i == s->held) ? 0.0 : s->inv_m[i];
}

static void build_chains(Sim *s);

int sim_rebuild(Sim *s)
{
    const Params *p = &s->p;
    const int C = s->cols, R = s->rows;
    const int per = p->bc == BC_PERIODIC;

    s->prestrain_eff = (p->bc == BC_CLAMPED || p->bc == BC_PERIODIC) ? p->prestrain : 0.0;
    s->nspr = 0;
    s->nrods = 0;

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
        if (s->spr[k].rod) continue;
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
            if (sp->rod) continue;
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

    build_chains(s);
    sim_update_damping(s);
    sim_compute_forces(s);

    /* bring the current state onto the rod constraints; if the rods cannot
     * accommodate it, fall back to the reference state (which satisfies
     * every rod) and report that */
    if (s->nrods && sim_project(s) > PROJ_TOL) {
        memset(s->ux, 0, (size_t)s->n * sizeof(double));
        memset(s->uy, 0, (size_t)s->n * sizeof(double));
        memset(s->vx, 0, (size_t)s->n * sizeof(double));
        memset(s->vy, 0, (size_t)s->n * sizeof(double));
        s->held = -1;
        sim_compute_forces(s);
        return 1;
    }
    return 0;
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
        if (sp->rod) continue;
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

/* ───────────────────────────────────────────────────────── rods ──
 *
 * On this lattice the rods of one family form straight chains (rows,
 * columns or diagonals).  Within a chain, consecutive rods share a node, so
 * the linearised constraint equations are tridiagonal and can be solved
 * exactly (Thomas algorithm).  The solver sweeps over all chains, family by
 * family (block Gauss-Seidel), re-linearising each sweep, until every rod
 * length is within ROD_TOL.  Rows and columns constrain different
 * displacement components, so a few sweeps suffice.
 *
 * Rod j of a chain joins node n_{j-1} (its a end) to n_j (its b end).
 */

static void build_chains(Sim *s)
{
    free(s->chain); free(s->chain_off);
    s->chain = NULL; s->chain_off = NULL; s->nchains = 0; s->max_chain = 0;
    if (!s->nrods) return;

    const int n = s->n;
    static const int dirs[4][2] = { { 1, 0 }, { 0, 1 }, { 1, 1 }, { -1, 1 } };
    int *from = malloc((size_t)n * sizeof *from), *into = malloc((size_t)n * sizeof *into);
    unsigned char *used = calloc((size_t)s->nspr, 1);
    s->chain = malloc((size_t)s->nrods * sizeof *s->chain);
    s->chain_off = malloc((size_t)(s->nrods + 1) * sizeof *s->chain_off);
    if (!from || !into || !used || !s->chain || !s->chain_off) { free(from); free(into); free(used); return; }

    int len = 0;
    for (int d = 0; d < 4; d++) {
        for (int i = 0; i < n; i++) from[i] = into[i] = -1;
        for (int k = 0; k < s->nspr; k++) {
            const Spring *sp = &s->spr[k];
            if (sp->rod && sp->rx == dirs[d][0] && sp->ry == dirs[d][1]) { from[sp->a] = k; into[sp->b] = k; }
        }
        /* open chains start at a node with an outgoing rod and no incoming one;
         * the second pass picks up closed rings (periodic) */
        for (int pass = 0; pass < 2; pass++)
            for (int i = 0; i < n; i++) {
                int k = from[i];
                if (k < 0 || used[k] || (pass == 0 && into[i] >= 0)) continue;
                s->chain_off[s->nchains++] = len;
                int start = i, node = i;
                do {
                    k = from[node];
                    if (k < 0 || used[k]) break;
                    used[k] = 1;
                    s->chain[len++] = k;
                    node = s->spr[k].b;
                } while (node != start);
                int cl = len - s->chain_off[s->nchains - 1];
                if (cl > s->max_chain) s->max_chain = cl;
            }
    }
    s->chain_off[s->nchains] = len;
    free(from); free(into); free(used);
    free(s->tri); free(s->vfac);
    s->tri  = malloc((size_t)(s->max_chain + 1) * 10 * sizeof *s->tri);
    s->vfac = malloc((size_t)s->nrods * 5 * sizeof *s->vfac);
}

/* Solve A lam = rhs for one tridiagonal chain system (Thomas algorithm).
 * lo, di, up, rhs have length m; lam is written into rhs. */
static void thomas(double *lo, double *di, double *up, double *rhs, int m)
{
    for (int j = 1; j < m; j++) {
        double f = lo[j] / di[j - 1];
        di[j]  -= f * up[j - 1];
        rhs[j] -= f * rhs[j - 1];
    }
    rhs[m - 1] /= di[m - 1];
    for (int j = m - 2; j >= 0; j--) rhs[j] = (rhs[j] - up[j] * rhs[j + 1]) / di[j];
}

static inline void bond(const Sim *s, const Spring *sp, const double *ux, const double *uy,
                        double *dx, double *dy)
{
    (void)s;
    *dx = sp->rx + ux[sp->b] - ux[sp->a];
    *dy = sp->ry + uy[sp->b] - uy[sp->a];
}

/*
 * Position solve.  Corrections move node a by +w_a lam e and node b by
 * -w_b lam e, with e the bond at the start of the step (SHAKE, dt > 0) or
 * the current bond (plain projection, dt = 0).  For rod j:
 *     sigma_j + 2 d_j . (dx_b - dx_a) = 0,  sigma_j = |d_j|^2 - L^2
 * couples lam_{j-1}, lam_j, lam_{j+1}.  With dt > 0 the same displacement
 * divided by dt goes into the velocities and the rod tension is recorded.
 */
static double solve_positions(Sim *s, double dt, int maxit, int *iters)
{
    const int use_old = dt > 0, M = s->max_chain + 1;
    double *lo = s->tri, *di = lo + M, *up = di + M, *rh = up + M;
    double *cx = rh + M, *cy = cx + M, *ex = cy + M, *ey = ex + M, *wA = ey + M, *wB = wA + M;
    double worst = 0;
    int it;
    if (use_old) for (int k = 0; k < s->nspr; k++) s->spr[k].f = 0.0;

    for (it = 0; it < maxit; it++) {
        worst = 0;
        for (int c = 0; c < s->nchains; c++) {
            const int *ch = &s->chain[s->chain_off[c]];
            const int m = s->chain_off[c + 1] - s->chain_off[c];
            double cworst = 0;
            for (int j = 0; j < m; j++) {
                const Spring *sp = &s->spr[ch[j]];
                bond(s, sp, s->ux, s->uy, &cx[j], &cy[j]);
                if (use_old) bond(s, sp, s->oux, s->ouy, &ex[j], &ey[j]);
                else { ex[j] = cx[j]; ey[j] = cy[j]; }
                wA[j] = wgt(s, sp->a); wB[j] = wgt(s, sp->b);
                double L2 = sp->L0 * sp->L0, sig = cx[j] * cx[j] + cy[j] * cy[j] - L2;
                double e = fabs(sig) / L2;
                if (!(e <= cworst)) cworst = isfinite(e) ? e : INFINITY;   /* NaN-safe */
                rh[j] = -sig;
            }
            if (!(cworst <= worst)) worst = cworst;
            if (cworst < ROD_TOL || !isfinite(cworst)) continue;
            for (int j = 0; j < m; j++) {
                if (wA[j] + wB[j] == 0) { lo[j] = up[j] = 0; di[j] = 1; rh[j] = 0; continue; }
                di[j] = -2.0 * (wA[j] + wB[j]) * (cx[j] * ex[j] + cy[j] * ey[j]);
                if (fabs(di[j]) < 1e-12) di[j] = -1e-12;
                lo[j] = j > 0     ? 2.0 * wA[j] * (cx[j] * ex[j - 1] + cy[j] * ey[j - 1]) : 0.0;
                up[j] = j < m - 1 ? 2.0 * wB[j] * (cx[j] * ex[j + 1] + cy[j] * ey[j + 1]) : 0.0;
            }
            thomas(lo, di, up, rh, m);
            for (int j = 0; j < m; j++) {
                Spring *sp = &s->spr[ch[j]];
                int a = sp->a, b = sp->b;
                double lam = rh[j], gx = lam * ex[j], gy = lam * ey[j];
                s->ux[a] += wA[j] * gx; s->uy[a] += wA[j] * gy;
                s->ux[b] -= wB[j] * gx; s->uy[b] -= wB[j] * gy;
                if (use_old) {
                    double idt = 1.0 / dt;
                    s->vx[a] += wA[j] * gx * idt; s->vy[a] += wA[j] * gy * idt;
                    s->vx[b] -= wB[j] * gx * idt; s->vy[b] -= wB[j] * gy * idt;
                    sp->f += 2.0 * lam * sqrt(ex[j] * ex[j] + ey[j] * ey[j]) * idt * idt;
                }
            }
        }
        if (worst < ROD_TOL || !isfinite(worst)) break;
    }
    if (iters) *iters = it + (it < maxit);
    return worst;
}

/* Velocity solve (RATTLE second half): make (v_b - v_a) . d_j = 0 for every
 * rod.  Positions are fixed here, so each chain's tridiagonal matrix is
 * constant: it is factored once and only the right-hand side changes from
 * sweep to sweep. */
static void solve_velocities(Sim *s)
{
    for (int c = 0; c < s->nchains; c++) {            /* factor every chain */
        const int *ch = &s->chain[s->chain_off[c]];
        const int m = s->chain_off[c + 1] - s->chain_off[c];
        double *F = &s->vfac[(size_t)5 * s->chain_off[c]];   /* per rod: dx dy lo' di' up */
        for (int j = 0; j < m; j++) {
            const Spring *sp = &s->spr[ch[j]];
            double *f = &F[5 * j];
            bond(s, sp, s->ux, s->uy, &f[0], &f[1]);
            double wa = wgt(s, sp->a), wb = wgt(s, sp->b);
            if (wa + wb == 0) { f[2] = 0; f[3] = 1; f[4] = 0; continue; }
            f[3] = -(wa + wb) * (f[0] * f[0] + f[1] * f[1]);
            f[2] = j > 0     ? wa * (f[0] * F[5 * (j - 1)] + f[1] * F[5 * (j - 1) + 1]) : 0.0;
            f[4] = 0.0;
        }
        for (int j = 0; j < m - 1; j++) {             /* up[j] needs d_{j+1} */
            const Spring *sp = &s->spr[ch[j]];
            double wb = wgt(s, sp->b), *f = &F[5 * j];
            if (wgt(s, sp->a) + wb == 0) continue;
            f[4] = wb * (f[0] * F[5 * (j + 1)] + f[1] * F[5 * (j + 1) + 1]);
        }
        for (int j = 1; j < m; j++) {                 /* forward elimination: lo' = l/di_prev */
            double *f = &F[5 * j], *p = &F[5 * (j - 1)];
            f[2] /= p[3];
            f[3] -= f[2] * p[4];
        }
    }

    const int M = s->max_chain + 1;
    double *rh = s->tri;
    (void)M;
    for (int it = 0; it < ROD_MAXIT; it++) {
        double worst = 0;
        for (int c = 0; c < s->nchains; c++) {
            const int *ch = &s->chain[s->chain_off[c]];
            const int m = s->chain_off[c + 1] - s->chain_off[c];
            const double *F = &s->vfac[(size_t)5 * s->chain_off[c]];
            double cworst = 0;
            for (int j = 0; j < m; j++) {
                const Spring *sp = &s->spr[ch[j]];
                const double *f = &F[5 * j];
                double cv = f[0] * (s->vx[sp->b] - s->vx[sp->a]) + f[1] * (s->vy[sp->b] - s->vy[sp->a]);
                double e = cv * cv / (f[0] * f[0] + f[1] * f[1]);
                if (e > cworst) cworst = e;
                rh[j] = -cv;
            }
            if (cworst > worst) worst = cworst;
            if (cworst < 1e-20) continue;
            for (int j = 1; j < m; j++) rh[j] -= F[5 * j + 2] * rh[j - 1];
            rh[m - 1] /= F[5 * (m - 1) + 3];
            for (int j = m - 2; j >= 0; j--) rh[j] = (rh[j] - F[5 * j + 4] * rh[j + 1]) / F[5 * j + 3];
            for (int j = 0; j < m; j++) {
                const Spring *sp = &s->spr[ch[j]];
                const double *f = &F[5 * j];
                int a = sp->a, b = sp->b;
                double wa = wgt(s, a), wb = wgt(s, b), gx = rh[j] * f[0], gy = rh[j] * f[1];
                s->vx[a] += wa * gx; s->vy[a] += wa * gy;
                s->vx[b] -= wb * gx; s->vy[b] -= wb * gy;
            }
        }
        if (worst < 1e-20) return;
    }
}

static void snapshot(Sim *s);
static void restore(Sim *s);

/* returns 0 on success; on failure (constraints cannot be met, e.g. a drag
 * the rods forbid) the caller restores the step-start state */
static int shake(Sim *s, double dt)
{
    s->rod_resid = solve_positions(s, dt, ROD_MAXIT, &s->rod_iters);
    /* an unconverged but small residual is fine: the next SHAKE targets the
     * exact length again, so the error does not accumulate */
    return (isfinite(s->rod_resid) && s->rod_resid < BLOCK_TOL) ? 0 : -1;
}

void sim_project_velocities(Sim *s)
{
    if (s->nrods) solve_velocities(s);
}

double sim_project(Sim *s)
{
    if (!s->nrods) return 0.0;
    double worst = solve_positions(s, 0.0, PROJ_MAXIT, NULL);
    if (!isfinite(worst)) return INFINITY;
    solve_velocities(s);
    sim_compute_forces(s);
    return worst;
}

void sim_step(Sim *s)
{
    const double dt = s->p.dt, h = 0.5 * dt;
    const int n = s->n, held = s->held;

    if (s->nrods) {
        snapshot(s);
        memcpy(s->oux, s->ux, (size_t)n * sizeof(double));
        memcpy(s->ouy, s->uy, (size_t)n * sizeof(double));
    }
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
    if (s->nrods && shake(s, dt)) {
        /* blocked: undo the step, stop the dragged node where it was */
        restore(s);
        for (int k = 0; k < s->nspr; k++) s->spr[k].f = 0.0;   /* failed solve: no valid forces */
        s->rod_blocked = 1;
        if (held >= 0) { s->held_ux = s->ux[held]; s->held_uy = s->uy[held]; }
        s->steps++;
        s->t = s->steps * dt;
        return;
    }
    s->rod_blocked = 0;
    sim_compute_forces(s);
    for (int i = 0; i < n; i++) {
        if (s->fixed[i] || i == held) continue;
        s->vx[i] += h * s->ax[i]; s->vy[i] += h * s->ay[i];
    }
    damp_half(s);
    if (held >= 0) { s->vx[held] = s->held_vx; s->vy[held] = s->held_vy; }
    if (s->nrods) solve_velocities(s);

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
        if (sp->rod) continue;
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
    L->tc = s->p.rod_card ? 0.0 : s->p.k * e;      /* T/|r| = k(|r| - L0)/|r| = k e */
    L->td = s->p.rod_diag ? 0.0 : s->p.kdiag * e;  /* rods carry no stored tension  */
    L->rod_card = s->p.rod_card; L->rod_diag = s->p.rod_diag;
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
    /* A kick should excite internal modes, not send the sheet drifting or
     * spinning off: remove net momentum where it is conserved (free,
     * periodic) and net angular momentum where that is conserved (free). */
    if (s->held < 0 && (s->p.bc == BC_FREE || s->p.bc == BC_PERIODIC)) {
        double M = 0, px = 0, py = 0, cx = 0, cy = 0;
        for (int i = 0; i < s->n; i++) {
            double x = i % s->cols + s->ux[i], y = i / s->cols + s->uy[i];
            M += s->m[i]; px += s->m[i] * s->vx[i]; py += s->m[i] * s->vy[i];
            cx += s->m[i] * x; cy += s->m[i] * y;
        }
        px /= M; py /= M; cx /= M; cy /= M;
        for (int i = 0; i < s->n; i++) { s->vx[i] -= px; s->vy[i] -= py; }
        if (s->p.bc == BC_FREE) {
            double Lz = 0, I = 0;
            for (int i = 0; i < s->n; i++) {
                double x = i % s->cols + s->ux[i] - cx, y = i / s->cols + s->uy[i] - cy;
                Lz += s->m[i] * (x * s->vy[i] - y * s->vx[i]);
                I  += s->m[i] * (x * x + y * y);
            }
            double W = Lz / I;
            for (int i = 0; i < s->n; i++) {
                double x = i % s->cols + s->ux[i] - cx, y = i / s->cols + s->uy[i] - cy;
                s->vx[i] += W * y; s->vy[i] -= W * x;
            }
        }
    }
    sim_project_velocities(s);
}

static void snapshot(Sim *s)
{
    size_t b = (size_t)s->n * sizeof(double);
    memcpy(s->sx, s->ux, b); memcpy(s->sy, s->uy, b);
    memcpy(s->svx, s->vx, b); memcpy(s->svy, s->vy, b);
}

static void restore(Sim *s)
{
    size_t b = (size_t)s->n * sizeof(double);
    memcpy(s->ux, s->sx, b); memcpy(s->uy, s->sy, b);
    memcpy(s->vx, s->svx, b); memcpy(s->vy, s->svy, b);
    sim_compute_forces(s);
}

int sim_displace(Sim *s, int i, double dx, double dy)
{
    if (i < 0 || i >= s->n || s->fixed[i]) return -1;
    snapshot(s);
    s->ux[i] += dx; s->uy[i] += dy;
    s->vx[i] = s->vy[i] = 0.0;
    if (s->nrods) {
        /* move the plucked node first, then let the rods drag its
         * neighbours along: hold it while projecting */
        int keep = s->held; s->held = i;
        double r = sim_project(s);
        s->held = keep;
        if (r > PROJ_TOL) { restore(s); return -1; }
    }
    sim_compute_forces(s);
    return 0;
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

int sim_excite_mode(Sim *s, int branch, double strain)
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
    if (!isfinite(w[branch])) {               /* branch frozen by rods */
        memset(s->ux, 0, (size_t)s->n * sizeof(double));
        memset(s->uy, 0, (size_t)s->n * sizeof(double));
        sim_compute_forces(s);
        return -1;
    }
    if (s->nrods && sim_project(s) > PROJ_TOL) {
        memset(s->ux, 0, (size_t)s->n * sizeof(double));
        memset(s->uy, 0, (size_t)s->n * sizeof(double));
        sim_compute_forces(s);
        return -1;
    }
    sim_compute_forces(s);
    return 0;
}
