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
                        &s->oux, &s->ouy, &s->sx, &s->sy, &s->svx, &s->svy,
                        &s->sax, &s->say, &s->noise };
    for (size_t i = 0; i < sizeof arrs / sizeof *arrs; i++)
        if (!(*arrs[i] = calloc(n, sizeof(double)))) return -1;
    if (!(s->fixed = calloc(n, 1))) return -1;
    s->spr_cap = s->n * 4 + 8;
    if (!(s->spr = calloc((size_t)s->spr_cap, sizeof *s->spr))) return -1;
    if (!(s->cell = calloc((size_t)s->n * 4, sizeof *s->cell))) return -1;
    if (!(s->exc = calloc((size_t)s->spr_cap, sizeof *s->exc))) return -1;
    if (!(s->tarea = calloc((size_t)s->n * 4, sizeof *s->tarea))) return -1;

    sim_set_masses(s);
    sim_rebuild(s);
    return 0;
}

void sim_free(Sim *s)
{
    free(s->ux); free(s->uy); free(s->vx); free(s->vy); free(s->ax); free(s->ay);
    free(s->m); free(s->inv_m); free(s->gamma); free(s->decay); free(s->anchor_k);
    free(s->oux); free(s->ouy); free(s->sx); free(s->sy); free(s->svx); free(s->svy);
    free(s->sax); free(s->say); free(s->noise); free(s->cell); free(s->exc); free(s->tarea);
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

static inline int kinematic(const Sim *s, int i) { return i == s->held && !s->drag_soft; }

static inline double wgt(const Sim *s, int i)
{
    return (s->fixed[i] || kinematic(s, i)) ? 0.0 : s->inv_m[i];
}

/* ── bond law: Hooke, or bounded with a wall at |r - L0| = D ───────────── */

/* energy and dV/dr of a (non-rod) bond at length r; INFINITY past the wall.
 * Bounded: V = k e^2 / (2g), g = 1 - (e/D)^2, so dV/de = k e / g^2. */
static inline void bond_law(const Sim *s, const Spring *sp, double r, double *V, double *dV)
{
    double e = r - sp->L0;
    if (s->p.bounded) {
        double D = s->p.bound_delta, x = e / D, g = 1.0 - x * x;
        if (g <= 0.0) { *V = INFINITY; *dV = e > 0 ? INFINITY : -INFINITY; return; }
        *V  = 0.5 * sp->k * e * e / g;
        *dV = sp->k * e / (g * g);
    } else {
        *V  = 0.5 * sp->k * e * e;
        *dV = sp->k * e;
    }
}

/* ── orientation barrier on the corner triangles of every unit cell ─────── */
/* b(A) = kappa (A_on - A)^2 / (A A_on) below A_on: C^1 there, 1/A at 0 */

static const double CORNER[4][2] = { { 0, 0 }, { 1, 0 }, { 1, 1 }, { 0, 1 } };

/* signed area of the triangle at corner k of cell c (positive at rest, 1/2) */
static inline double tri_area(const Sim *s, const int *cn, int k, double g[3][2])
{
    int k1 = (k + 1) & 3, k3 = (k + 3) & 3;
    int i0 = cn[k], i1 = cn[k1], i2 = cn[k3];
    double ax = CORNER[k1][0] - CORNER[k][0] + s->ux[i1] - s->ux[i0];
    double ay = CORNER[k1][1] - CORNER[k][1] + s->uy[i1] - s->uy[i0];
    double bx = CORNER[k3][0] - CORNER[k][0] + s->ux[i2] - s->ux[i0];
    double by = CORNER[k3][1] - CORNER[k][1] + s->uy[i2] - s->uy[i0];
    if (g) {                                  /* dA/dp for p0, p1, p2 */
        g[1][0] =  0.5 * by;  g[1][1] = -0.5 * bx;
        g[2][0] = -0.5 * ay;  g[2][1] =  0.5 * ax;
        g[0][0] = -g[1][0] - g[2][0];  g[0][1] = -g[1][1] - g[2][1];
    }
    return 0.5 * (ax * by - bx * ay);
}

static inline void barrier(const Sim *s, double A, double *b, double *db)
{
    const double Aon = 0.5 * s->p.orient_floor, kap = s->p.orient_kappa;
    if (A >= Aon) { *b = 0; *db = 0; return; }
    if (A <= 0)   { *b = INFINITY; *db = -INFINITY; return; }
    double d = Aon - A;
    *b  =  kap * d * d / (A * Aon);
    *db = -kap * (Aon * Aon - A * A) / (Aon * A * A);
}

static void build_cells(Sim *s)
{
    const int C = s->cols, R = s->rows, per = s->p.bc == BC_PERIODIC;
    s->ncells = 0;
    for (int r = 0; r < (per ? R : R - 1); r++)
        for (int c = 0; c < (per ? C : C - 1); c++) {
            int c1 = (c + 1) % C, r1 = (r + 1) % R;
            int *cn = &s->cell[4 * s->ncells++];
            cn[0] = sim_idx(s, c, r);  cn[1] = sim_idx(s, c1, r);
            cn[2] = sim_idx(s, c1, r1); cn[3] = sim_idx(s, c, r1);
        }
}

int sim_valid(const Sim *s)
{
    if (s->p.bounded) {
        const double D = s->p.bound_delta;
        for (int k = 0; k < s->nspr; k++) {
            const Spring *sp = &s->spr[k];
            if (sp->rod) continue;
            double dx = sp->rx + s->ux[sp->b] - s->ux[sp->a], dy = sp->ry + s->uy[sp->b] - s->uy[sp->a];
            double e = sqrt(dx * dx + dy * dy) - sp->L0;
            if (!(fabs(e) < D * (1.0 - 1e-9))) return 0;
        }
    }
    if (s->p.orient)
        for (int c = 0; c < s->ncells; c++)
            for (int k = 0; k < 4; k++)
                if (!(tri_area(s, &s->cell[4 * c], k, NULL) > 1e-12)) return 0;
    return 1;
}

/* Was the substep that moved the lattice from (oux, ouy) to (ux, uy) fine
 * enough near the walls?  Every bond may close at most half of its gap to
 * the bound, and every triangle below A_on may lose at most half its area.
 * This also implies validity, and it is what keeps Verlet accurate (and
 * energy-conserving) when a bond comes close to its wall. */
static int step_resolved(const Sim *s)
{
    if (s->p.bounded) {
        const double D = s->p.bound_delta;
        for (int k = 0; k < s->nspr; k++) {
            const Spring *sp = &s->spr[k];
            if (sp->rod) continue;
            double x1 = sp->rx + s->ux[sp->b] - s->ux[sp->a], y1 = sp->ry + s->uy[sp->b] - s->uy[sp->a];
            double e0 = s->exc[k], e1 = sqrt(x1 * x1 + y1 * y1) - sp->L0;
            if (!(fabs(e1 - e0) <= 0.5 * (D - fabs(e0)))) return 0;
        }
    }
    if (s->p.orient) {
        const double Aon = 0.5 * s->p.orient_floor;
        for (int c = 0; c < s->ncells; c++)
            for (int k = 0; k < 4; k++) {
                double A1 = tri_area(s, &s->cell[4 * c], k, NULL);
                if (A1 < Aon && !(A1 > 0.5 * s->tarea[4 * c + k])) return 0;
            }
    }
    return 1;
}

/* Highest local frequency set by the walls: sqrt(V'' (1/m_a + 1/m_b)) for
 * bounded bonds, sqrt(b'' sum_j |dA/dp_j|^2 / m_j) for triangles below A_on,
 * computed during the force pass.  Chooses the number of substeps before a
 * step is attempted. */
static double wall_omega(const Sim *s) { return sqrt(s->wall_w2); }   /* cached by the force pass */

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
        double len = sqrt(s->spr[k].rx * s->spr[k].rx + s->spr[k].ry * s->spr[k].ry), V, dV;
        bond_law(s, &s->spr[k], len, &V, &dV);
        pe += V;
    }
    s->pe_ref = pe;
    build_cells(s);

    /* Gershgorin bound on omega_max^2 */
    double *row = calloc((size_t)s->n, sizeof *row);
    if (row) {
        for (int k = 0; k < s->nspr; k++) {
            const Spring *sp = &s->spr[k];
            if (sp->rod) continue;
            double len = sqrt(sp->rx * sp->rx + sp->ry * sp->ry);
            double t = sp->k * (len - sp->L0) / len, keff = sp->k;
            if (s->p.bounded) {
                double x = (len - sp->L0) / s->p.bound_delta, g = 1 - x * x;
                keff = sp->k * (1 + 3 * x * x) / (g * g * g);
                t /= g * g;
            }
            double w = 2.0 * (keff + fabs(t));
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
    if ((s->nrods && sim_project(s) > PROJ_TOL) || !sim_valid(s)) {
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
            if (p->thermo) g += p->therm_gamma;
            s->gamma[i] = g;
            s->decay[i] = exp(-g * p->dt * 0.5 * s->inv_m[i]);
            s->noise[i] = p->thermo ? sqrt(p->temp * s->inv_m[i] * (1.0 - s->decay[i] * s->decay[i])) : 0.0;
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

    double wmax2 = 0;
    for (int i = 0; i < n; i++) {
        ax[i] = -s->anchor_k[i] * ux[i];
        ay[i] = -s->anchor_k[i] * uy[i];
    }
    if (s->held >= 0 && s->drag_soft) {       /* mouse spring */
        int h = s->held;
        double kd = 10.0 * s->p.k, cd = 2.0 * sqrt(kd * s->m[h]);
        double dx = s->held_ux - ux[h], dy = s->held_uy - uy[h], d = sqrt(dx * dx + dy * dy);
        if (d > 1.0) { dx /= d; dy /= d; }
        ax[h] += kd * dx + cd * (s->held_vx - s->vx[h]);
        ay[h] += kd * dy + cd * (s->held_vy - s->vy[h]);
    }
    for (int k = 0; k < s->nspr; k++) {
        const Spring *sp = &s->spr[k];
        if (sp->rod) continue;
        int a = sp->a, b = sp->b;
        double dx = sp->rx + ux[b] - ux[a];
        double dy = sp->ry + uy[b] - uy[a];
        double dist = sqrt(dx * dx + dy * dy);
        if (dist < 1e-12) continue;
        double V, dV;
        bond_law(s, sp, dist, &V, &dV);
        if (!isfinite(dV)) continue;             /* past a wall: the step is rejected anyway */
        if (s->p.bounded) {
            double e = dist - sp->L0, x = e / s->p.bound_delta, g = 1 - x * x;
            double w2 = sp->k * (1 + 3 * x * x) / (g * g * g) * (s->inv_m[a] + s->inv_m[b]);
            s->exc[k] = e;
            if (w2 > wmax2) wmax2 = w2;
        }
        double f = dV / dist;
        double fx = f * dx, fy = f * dy;
        ax[a] += fx; ay[a] += fy;
        ax[b] -= fx; ay[b] -= fy;
    }
    if (s->p.orient)
        for (int c = 0; c < s->ncells; c++) {
            const int *cn = &s->cell[4 * c];
            for (int k = 0; k < 4; k++) {
                double g[3][2], b, db;
                double A = tri_area(s, cn, k, g);
                s->tarea[4 * c + k] = A;
                barrier(s, A, &b, &db);
                if (db == 0.0 || !isfinite(db)) continue;
                int id[3] = { cn[k], cn[(k + 1) & 3], cn[(k + 3) & 3] };
                double sum = 0;
                for (int j = 0; j < 3; j++) {
                    ax[id[j]] -= db * g[j][0]; ay[id[j]] -= db * g[j][1];
                    sum += (g[j][0] * g[j][0] + g[j][1] * g[j][1]) * s->inv_m[id[j]];
                }
                double w2 = 2.0 * s->p.orient_kappa * 0.5 * s->p.orient_floor / (A * A * A) * sum;
                if (w2 > wmax2) wmax2 = w2;
            }
        }
    s->wall_w2 = wmax2;
    for (int i = 0; i < n; i++) { ax[i] *= s->inv_m[i]; ay[i] *= s->inv_m[i]; }
}

/* O half-step: exact drag over dt/2, plus thermostat noise.  For the
 * nominal dt the factors are precomputed; substeps compute their own. */
static void damp_half(Sim *s, double dt)
{
    const int nominal = dt == s->p.dt, thermo = s->p.thermo;
    for (int i = 0; i < s->n; i++) {
        if (s->fixed[i]) continue;
        double c = nominal ? s->decay[i] : exp(-s->gamma[i] * dt * 0.5 * s->inv_m[i]);
        s->vx[i] *= c; s->vy[i] *= c;
        if (thermo) {
            /* unit-variance uniform noise, sqrt(3) U[-1,1]: each step adds a
             * tiny kick (1 - c ~ 1e-3), so over the relaxation time the sum is
             * Gaussian by the central limit theorem; equipartition is tested */
            double sg = (nominal ? s->noise[i] : sqrt(s->p.temp * s->inv_m[i] * (1.0 - c * c))) * 1.7320508075688772;
            s->vx[i] += sg * rng_sym(&s->rng);
            s->vy[i] += sg * rng_sym(&s->rng);
        }
    }
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

#define SUB_MAXDEPTH 10         /* up to 1024 substeps per step (only near walls) */

/* one integration step of length dt; returns -1 if it would cross a wall or
 * SHAKE cannot close (the caller restores the step-start state) */
static int substep(Sim *s, double dt)
{
    const double h = 0.5 * dt;
    const int n = s->n, held = s->held;

    if (s->nrods || s->p.bounded || s->p.orient) {
        memcpy(s->oux, s->ux, (size_t)n * sizeof(double));
        memcpy(s->ouy, s->uy, (size_t)n * sizeof(double));
    }
    damp_half(s, dt);
    for (int i = 0; i < n; i++) {
        if (s->fixed[i] || kinematic(s, i)) continue;
        s->vx[i] += h * s->ax[i];  s->vy[i] += h * s->ay[i];
        s->ux[i] += dt * s->vx[i]; s->uy[i] += dt * s->vy[i];
    }
    if (held >= 0 && !s->drag_soft) {
        s->ux[held] = s->held_ux; s->uy[held] = s->held_uy;
        s->vx[held] = s->held_vx; s->vy[held] = s->held_vy;
    }
    if (s->nrods && shake(s, dt)) return -1;
    if ((s->p.bounded || s->p.orient) && !step_resolved(s)) return -1;
    sim_compute_forces(s);
    for (int i = 0; i < n; i++) {
        if (s->fixed[i] || kinematic(s, i)) continue;
        s->vx[i] += h * s->ax[i]; s->vy[i] += h * s->ay[i];
    }
    damp_half(s, dt);
    if (held >= 0 && !s->drag_soft) { s->vx[held] = s->held_vx; s->vy[held] = s->held_vy; }
    if (s->nrods) solve_velocities(s);
    return 0;
}

static void save_full(Sim *s)
{
    size_t b = (size_t)s->n * sizeof(double);
    memcpy(s->sx, s->ux, b);  memcpy(s->sy, s->uy, b);
    memcpy(s->svx, s->vx, b); memcpy(s->svy, s->vy, b);
    memcpy(s->sax, s->ax, b); memcpy(s->say, s->ay, b);
}

static void restore_full(Sim *s)
{
    size_t b = (size_t)s->n * sizeof(double);
    memcpy(s->ux, s->sx, b);  memcpy(s->uy, s->sy, b);
    memcpy(s->vx, s->svx, b); memcpy(s->vy, s->svy, b);
    memcpy(s->ax, s->sax, b); memcpy(s->ay, s->say, b);
}

void sim_step(Sim *s)
{
    const double dt = s->p.dt;
    const int guarded = s->nrods || s->p.bounded || s->p.orient;
    int depth = 0, ok = 1;

    int depth0 = 0;
    if (guarded) {
        save_full(s);
        /* resolve the stiffest wall: omega_loc * dt_sub <= 0.1 */
        if (s->p.bounded || s->p.orient) {
            double need = wall_omega(s) * dt / 0.1;
            while (depth0 < SUB_MAXDEPTH && (double)(1 << depth0) < need) depth0++;
        }
    }
    for (depth = depth0; depth <= (guarded ? SUB_MAXDEPTH : 0); depth++) {
        int nsub = 1 << depth;
        ok = 1;
        for (int k = 0; k < nsub && ok; k++) if (substep(s, dt / nsub)) ok = 0;
        if (ok) break;
        restore_full(s);
        if (s->p.bounded || s->p.orient) sim_compute_forces(s);   /* refresh the wall caches */
    }
    s->sub_depth = ok ? depth : SUB_MAXDEPTH;
    s->rod_blocked = !ok;
    if (!ok) {
        /* even 1024 substeps cannot resolve the wall.  A drag into an impossible
         * configuration is refused: keep the step-start state and stop the
         * dragged node where it is.  Free motion must not freeze: reverse
         * every velocity (an exact time reversal, energy unchanged), so the
         * lattice retraces its path away from the wall. */
        for (int k = 0; k < s->nspr; k++) s->spr[k].f = 0.0;
        if (s->held >= 0 && !s->drag_soft) { s->held_ux = s->ux[s->held]; s->held_uy = s->uy[s->held]; }
        else {
            for (int i = 0; i < s->n; i++) { s->vx[i] = -s->vx[i]; s->vy[i] = -s->vy[i]; }
            s->reversals++;
            s->rod_blocked = 0;
        }
    }
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
        double V, dV;
        bond_law(s, sp, sqrt(dx * dx + dy * dy), &V, &dV);
        v += V;
    }
    if (s->p.orient)
        for (int c = 0; c < s->ncells; c++)
            for (int q = 0; q < 4; q++) {
                double b, db;
                barrier(s, tri_area(s, &s->cell[4 * c], q, NULL), &b, &db);
                v += b;
            }
    *ke = k; *pe = v - s->pe_ref;
}

double sim_lindemann(const Sim *s)
{
    double sum = 0; int n = 0;
    for (int j = 0; j < s->nspr; j++) {
        const Spring *sp = &s->spr[j];
        if (sp->rx != 0 && sp->ry != 0) continue;          /* nearest neighbours only */
        double dx = s->ux[sp->b] - s->ux[sp->a], dy = s->uy[sp->b] - s->uy[sp->a];
        sum += dx * dx + dy * dy; n++;
    }
    return n ? sqrt(sum / n) : 0.0;
}

void sim_lattice(const Sim *s, Lattice *L)
{
    double e = s->prestrain_eff;
    double kc = s->p.k, kd = s->p.kdiag, tc = s->p.k * e, td = s->p.kdiag * e;
    if (s->p.bounded) {
        /* linearise V about the reference length: k_eff = V'', T = V' */
        double D = s->p.bound_delta;
        double xc = e / D, gc = 1 - xc * xc;                   /* |r| = 1      */
        double xd = sqrt(2.0) * e / D, gd = 1 - xd * xd;       /* |r| = sqrt 2 */
        kc *= (1 + 3 * xc * xc) / (gc * gc * gc);  tc /= gc * gc;
        kd *= (1 + 3 * xd * xd) / (gd * gd * gd);  td /= gd * gd;
    }
    L->k = kc; L->kd = kd; L->diag = s->p.diagonal;
    L->tc = s->p.rod_card ? 0.0 : tc;      /* T/|r| */
    L->td = s->p.rod_diag ? 0.0 : td;      /* rods carry no stored tension */
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
    if (s->p.bounded) {                          /* keep a pluck inside the wall */
        double d = sqrt(dx * dx + dy * dy), cap = 0.8 * s->p.bound_delta;
        if (d > cap) { dx *= cap / d; dy *= cap / d; }
    }
    snapshot(s);
    s->ux[i] += dx; s->uy[i] += dy;
    s->vx[i] = s->vy[i] = 0.0;
    if (s->nrods) {
        /* move the plucked node first, then let the rods drag its
         * neighbours along: hold it while projecting */
        int keep = s->held, soft = s->drag_soft;
        s->held = i; s->drag_soft = 0;
        double r = sim_project(s);
        s->held = keep; s->drag_soft = soft;
        if (r > PROJ_TOL) { restore(s); return -1; }
    }
    if (!sim_valid(s)) { restore(s); return -1; }
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
    if ((s->nrods && sim_project(s) > PROJ_TOL) || !sim_valid(s)) {
        memset(s->ux, 0, (size_t)s->n * sizeof(double));
        memset(s->uy, 0, (size_t)s->n * sizeof(double));
        sim_compute_forces(s);
        return -1;
    }
    sim_compute_forces(s);
    return 0;
}
