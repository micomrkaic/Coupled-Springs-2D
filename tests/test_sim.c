/*
 * test_sim.c — headless checks of the physics and the spectral analysis
 *
 *   make test
 */
#include "../src/modes.h"
#include "../src/params.h"
#include "../src/sim.h"
#include "../src/spectrum.h"
#include "../src/util.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static int failures = 0;

#define CHECK(cond, ...) do {                                   \
    if (cond) { printf("  ok    "); printf(__VA_ARGS__); }      \
    else      { printf("  FAIL  "); printf(__VA_ARGS__); failures++; } \
    printf("\n"); } while (0)

static Params base(int n, Boundary bc, int diag)
{
    Params p; params_default(&p);
    p.cols = p.rows = n; p.bc = bc; p.diagonal = diag;
    p.kdiag = p.k; p.absorb_gamma = sqrt(p.k * p.mass);
    p.damping = 0.0; p.seed = 12345;
    return p;
}

/* 1. D(q) from modes.c equals the linearised force of sim.c */
static void test_dynmat(void)
{
    printf("dynamical matrix vs. linearised forces\n");
    const int N = 16;
    for (int variant = 0; variant < 4; variant++) {
        Params p = base(N, BC_PERIODIC, variant > 0);
        p.prestrain = variant >= 2 ? 0.07 : 0.0;
        p.bounded = variant == 3;                /* bounded: k_eff = V'', T = V' */
        Sim s; sim_init(&s, &p);
        Lattice L; sim_lattice(&s, &L);
        double worst = 0;
        for (int mx = 0; mx <= N / 2; mx += 3)
            for (int my = -N / 2 + 1; my <= N / 2; my += 5) {
                double qx = 2 * PI * mx / N, qy = 2 * PI * my / N;
                for (int pol = 0; pol < 2; pol++) {
                    const double eps = 1e-6;
                    for (int r = 0; r < N; r++)
                        for (int c = 0; c < N; c++) {
                            int i = sim_idx(&s, c, r);
                            double ph = cos(qx * c + qy * r);
                            s.ux[i] = pol == 0 ? eps * ph : 0;
                            s.uy[i] = pol == 1 ? eps * ph : 0;
                        }
                    sim_compute_forces(&s);
                    double D[3]; dynmat(&L, qx, qy, D);
                    /* a = -D u at a site with cos = 1 (c = r = 0) */
                    double ex = pol == 0 ? -D[0] * eps : -D[2] * eps;
                    double ey = pol == 0 ? -D[2] * eps : -D[1] * eps;
                    double err = fabs(s.ax[0] - ex) + fabs(s.ay[0] - ey);
                    double scale = (fabs(D[0]) + fabs(D[1]) + fabs(D[2]) + 1) * eps;
                    if (err / scale > worst) worst = err / scale;
                }
            }
        CHECK(worst < 1e-4, "variant %d (diag=%d, prestrain=%.2f, bounded=%d): max rel err %.2e",
              variant, p.diagonal, p.prestrain, p.bounded, worst);
        sim_free(&s);
    }
}

/* 2. energy is conserved (bounded error) with gamma = 0 */
static void test_energy(void)
{
    printf("energy conservation, gamma = 0, 20 s of sim time\n");
    Boundary bcs[] = { BC_FREE, BC_PINNED, BC_CLAMPED, BC_PERIODIC };
    for (int k = 0; k < 4; k++) {
        Params p = base(32, bcs[k], 1);
        p.prestrain = 0.05;
        Sim s; sim_init(&s, &p);
        sim_kick_random(&s, 2.0, 2.0);
        double ke, pe; sim_energy(&s, &ke, &pe);
        double E0 = ke + pe, worst = 0, last_q = 0, first_q = 0;
        int steps = 40000;
        for (int i = 0; i < steps; i++) {
            sim_step(&s);
            if (i % 100 == 0) {
                sim_energy(&s, &ke, &pe);
                double d = fabs(ke + pe - E0) / E0;
                if (d > worst) worst = d;
                if (i < steps / 4 && d > first_q) first_q = d;
                if (i > 3 * steps / 4 && d > last_q) last_q = d;
            }
        }
        CHECK(worst < 2e-4 && last_q < 4 * first_q + 1e-7,
              "%-8s max |dE|/E = %.2e (first quarter %.2e, last quarter %.2e)",
              bc_name(bcs[k]), worst, first_q, last_q);
        sim_free(&s);
    }
}

/* 3. operator-split drag is exact for a rigid translation */
static void test_damping(void)
{
    printf("exact drag on a rigid translation\n");
    Params p = base(8, BC_FREE, 1);
    p.damping = 0.7;
    Sim s; sim_init(&s, &p);
    for (int i = 0; i < s.n; i++) { s.vx[i] = 1.5; s.vy[i] = -0.5; }
    for (int i = 0; i < 4000; i++) sim_step(&s);
    double expect = 1.5 * exp(-p.damping * s.t / p.mass);
    double err = fabs(s.vx[0] - expect) / expect;
    CHECK(err < 1e-12, "v(t)/v0 vs exp(-gamma t/m) after %.1f s: rel err %.1e", s.t, err);
    sim_free(&s);
}

/* 4. peaks of the space-time spectrum sit on the analytic branches */
/* exact = 1: the basis diagonalises the linear problem, so every peak must
 * sit within 1.5 bins.  exact = 0 (free edges + diagonals): the DCT is only
 * approximate and free surfaces add Rayleigh-like edge waves, so require
 * the median relative error to be under 5%. */
static int cmp_dbl(const void *a, const void *b)
{ double x = *(const double *)a, y = *(const double *)b; return (x > y) - (x < y); }

static void run_dispersion_r(Boundary bc, int diag, double prestrain, int rods, int exact, double wscale,
                             double thr, const char *label)
{
    Params p = base(32, bc, diag);
    p.prestrain = prestrain;
    p.rod_card = rods & 1; p.rod_diag = (rods >> 1) & 1;
    p.nt = 512;
    Sim s; sim_init(&s, &p);
    Spectrum sp; spectrum_init(&sp, &s);
    /* small random displacements, zero velocity: a velocity kick would set
     * the floppy modes of the cardinal-only free lattice drifting, and
     * their second-order stretching shows up as low-frequency power */
    for (int i = 0; i < s.n; i++)
        if (!s.fixed[i]) { s.ux[i] = 1e-4 * rng_sym(&s.rng); s.uy[i] = 1e-4 * rng_sym(&s.rng); }
    sim_project(&s);
    sim_compute_forces(&s);
    while (sp.nsamp < sp.nt) { sim_step(&s); spectrum_after_step(&sp, &s); }
    int w, h; spectrum_disp_image(&sp, POL_L, &w, &h);   /* forces compute */

    double bin = sp.omega_nyq / (sp.nt / 2);
    int checked = 0, bad = 0; double worst = 0, rel[2048];
    for (int b = 0; b < 2; b++)
        for (int ix = 1; ix < sp.nlq; ix++) {
            double wa = (b == 0 ? sp.wL[ix] : sp.wT[ix]) * wscale;
            if (!isfinite(wa) || wa < 3 * bin) continue;  /* frozen, or too close to DC */
            int best = 0; double pmax = -1;
            for (int iw = 1; iw <= sp.nt / 2; iw++) {
                double v = sp.dpow[b][(size_t)iw * sp.nlq + ix];
                if (v > pmax) { pmax = v; best = iw; }
            }
            double err = fabs(best * bin - wa) / bin;
            if (err > worst) worst = err;
            if (err > 1.5) bad++;
            rel[checked++] = fabs(best * bin - wa) / wa;
        }
    qsort(rel, (size_t)checked, sizeof *rel, cmp_dbl);
    double median = checked ? rel[checked / 2] : 1.0;
    if (exact)
        CHECK(bad == 0 && checked > 0, "%-38s %d peaks, worst |dw| = %.2f bins",
              label, checked, worst);
    else
        CHECK(median < thr && checked > 0, "%-38s %d peaks, median |dw|/w = %.2f%% (limit %.1f%%)",
              label, checked, 100 * median, 100 * thr);
    spectrum_free(&sp); sim_free(&s);
}

static void run_dispersion(Boundary bc, int diag, double prestrain, int exact, const char *label)
{
    run_dispersion_r(bc, diag, prestrain, 0, exact, 1.0, 0.05, label);
}

static void test_dispersion(void)
{
    printf("omega-q peaks vs analytic branches (qy = lowest line)\n");
    run_dispersion(BC_PERIODIC, 0, 0.00, 1, "periodic, cardinal");
    run_dispersion(BC_PERIODIC, 1, 0.00, 1, "periodic, + diagonal");
    run_dispersion(BC_PERIODIC, 1, 0.08, 1, "periodic, + diagonal, prestrain 8%");
    run_dispersion(BC_FREE,     0, 0.00, 1, "free (DCT), cardinal");
    run_dispersion(BC_FREE,     1, 0.00, 0, "free (DCT), + diagonal");
    run_dispersion(BC_CLAMPED,  0, 0.08, 1, "clamped (DST), cardinal, prestrain 8%");
    /* free + cardinal rods: each column is rigid and includes the two edge
     * nodes, which have half the diagonal bonds, so the column stiffness is
     * (R-1) k_d rather than R k_d: omega = sqrt((R-1)/R) x bulk value */
    run_dispersion_r(BC_FREE,   1, 0.00, 1, 0, sqrt(31.0 / 32.0), 0.006, "free, card. RODS, w x sqrt((R-1)/R)");
    run_dispersion_r(BC_FREE,   1, 0.00, 1, 0, 1.0, 0.05, "  same, bulk formula (expect ~1.6%)");
    run_dispersion_r(BC_PINNED, 1, 0.00, 1, 0, 1.0, 0.05, "pinned, cardinal RODS + diag springs");
}

/* rods: constraints hold, energy is conserved, blocked moves are refused */
static double max_rod_err(const Sim *s)
{
    double worst = 0;
    for (int k = 0; k < s->nspr; k++) {
        const Spring *sp = &s->spr[k];
        if (!sp->rod) continue;
        double dx = sp->rx + s->ux[sp->b] - s->ux[sp->a], dy = sp->ry + s->uy[sp->b] - s->uy[sp->a];
        double e = fabs(sqrt(dx * dx + dy * dy) - sp->L0) / sp->L0;
        if (e > worst) worst = e;
    }
    return worst;
}

static void test_rods(void)
{
    printf("rigid rods (RATTLE)\n");
    struct { int card, diag, bc; const char *name; } cases[] = {
        { 1, 0, BC_FREE,   "free, cardinal rods + diag springs" },
        { 0, 1, BC_FREE,   "free, diagonal rods + card springs" },
        { 1, 0, BC_PINNED, "pinned, cardinal rods + diag springs" },
    };
    for (int c = 0; c < 3; c++) {
        Params p = base(32, cases[c].bc, 1);
        p.rod_card = cases[c].card; p.rod_diag = cases[c].diag;
        Sim s; sim_init(&s, &p);
        sim_displace(&s, sim_idx(&s, 16, 16), 0.0, -1.0);
        sim_kick_random(&s, 1.0, 1.0);
        double ke, pe; sim_energy(&s, &ke, &pe);
        double E0 = ke + pe, worst = 0, rerr = 0; int maxit = 0;
        for (int i = 0; i < 40000; i++) {
            sim_step(&s);
            if (s.rod_iters > maxit) maxit = s.rod_iters;
            if (i % 200 == 0) {
                sim_energy(&s, &ke, &pe);
                double d = fabs(ke + pe - E0) / E0; if (d > worst) worst = d;
                double r = max_rod_err(&s); if (r > rerr) rerr = r;
            }
        }
        CHECK(worst < 2e-4 && rerr < 1e-8, "%-38s |dE|/E %.1e, rod length err %.1e, max SHAKE sweeps %d",
              cases[c].name, worst, rerr, maxit);
        sim_free(&s);
    }
    /* clamped edges + cardinal rods: every row is a taut chain -> locked */
    Params p = base(16, BC_CLAMPED, 1);
    p.rod_card = 1;
    Sim s; sim_init(&s, &p);
    int r = sim_displace(&s, sim_idx(&s, 8, 8), 0.0, -0.5);
    double u = 0; for (int i = 0; i < s.n; i++) u += fabs(s.ux[i]) + fabs(s.uy[i]);
    CHECK(r == -1 && u == 0.0, "clamped + cardinal rods: pluck refused (%d), state untouched", r);
    /* a drag the rods forbid: the step is undone, nothing moves, no NaN */
    s.held = sim_idx(&s, 8, 8); s.held_ux = 0.0; s.held_uy = -0.3; s.held_vx = s.held_vy = 0.0;
    int blocked = 0;
    for (int i = 0; i < 100; i++) { sim_step(&s); blocked |= s.rod_blocked; }
    u = 0; for (int i = 0; i < s.n; i++) u += fabs(s.ux[i]) + fabs(s.uy[i]);
    CHECK(blocked && u == 0.0, "clamped + cardinal rods: forbidden drag blocked, state finite");
    sim_free(&s);
}

/* 5. a small clamped T mode stays in its mode */
static void test_mode_purity(void)
{
    printf("clamped string mode stays pure at small amplitude\n");
    Params p = base(32, BC_CLAMPED, 0);
    p.prestrain = 0.05;
    Sim s; sim_init(&s, &p);
    Spectrum sp; spectrum_init(&sp, &s);
    sim_excite_mode(&s, 1, 1e-4);
    for (int i = 0; i < 20000; i++) { sim_step(&s); spectrum_after_step(&sp, &s); }
    const double *v = &sp.modes.v[(size_t)(sp.modes.len - 1) * sp.modes.nser];
    double tot = 0; for (int k = 0; k < 2 * SPEC_KMODES; k++) tot += v[k];
    double frac = v[SPEC_KMODES + 0] / tot;
    CHECK(frac > 0.999, "energy fraction in T w=1 after %.0f s: %.6f", s.t, frac);
    spectrum_free(&sp); sim_free(&s);
}

/* 5b. FPUT recurrence: periodic lattice, diagonals, lowest L mode uniform
 * in y.  By symmetry every row moves identically, so this is a 1-D chain
 * whose diagonal bonds supply the anharmonicity. */
static void test_fput(void)
{
    printf("FPUT recurrence (32x32 periodic, diagonals, L mode, strain 0.15)\n");
    Params p = base(32, BC_PERIODIC, 1);
    Sim s; sim_init(&s, &p);
    Spectrum sp; spectrum_init(&sp, &s);
    sim_excite_mode(&s, 0, 0.15);
    spectrum_clear_history(&sp, &s);
    double ke, pe; sim_energy(&s, &ke, &pe);
    double E0 = ke + pe, fmin = 1, tmin = 0, fback = 0, tback = 0;
    while (s.t < 70.0) {
        sim_step(&s); spectrum_after_step(&sp, &s);
        if (s.steps % 200) continue;
        const double *v = &sp.modes.v[(size_t)(sp.modes.len - 1) * sp.modes.nser];
        double f = v[0] / E0;
        if (f < fmin) { fmin = f; tmin = s.t; }
        if (fmin < 0.4 && f > fback) { fback = f; tback = s.t; }
    }
    CHECK(fmin < 0.4 && fback > 0.9, "w=1 share falls to %.2f (t=%.0f s), returns to %.2f (t=%.0f s)",
          fmin, tmin, fback, tback);
    spectrum_free(&sp); sim_free(&s);
}

/* 6. light defect has a localised mode above the band */
static void test_defect(void)
{
    printf("light mass defect: local mode above the band\n");
    Params p = base(32, BC_FREE, 1);
    p.disorder = DIS_LIGHT; p.defect_light = 0.2;
    Sim s; sim_init(&s, &p);
    sim_set_masses(&s); sim_rebuild(&s);
    int d = sim_defect_index(&s);
    s.vx[d] = 1e-3;
    const int NS = 8192, stride = 4;
    double *x = malloc(NS * sizeof *x);
    for (int i = 0; i < NS; i++) { for (int j = 0; j < stride; j++) sim_step(&s); x[i] = s.ux[d]; }
    double dts = stride * p.dt, best_w = 0, best_p = -1;
    for (int k = 1; k < NS / 2; k += 1) {
        double wk = 2 * PI * k / (NS * dts), re = 0, im = 0;
        if (wk > 200) break;
        for (int i = 0; i < NS; i++) { re += x[i] * cos(wk * i * dts); im += x[i] * sin(wk * i * dts); }
        double pw = re * re + im * im;
        if (pw > best_p) { best_p = pw; best_w = wk; }
    }
    Lattice L; sim_lattice(&s, &L);
    L.m = p.mass;
    double top = band_top(&L);
    CHECK(best_w > 1.05 * top, "dominant defect frequency %.1f rad/s, band top %.1f rad/s", best_w, top);
    free(x); sim_free(&s);
}

/* 7. absorbing layer removes outgoing waves */
static double absorb_residual(int absorb, double gamma)
{
    Params p = base(64, BC_FREE, 1);
    p.absorb = absorb; p.absorb_width = 10; p.absorb_gamma = gamma;
    Sim s; sim_init(&s, &p);
    /* smooth pulse (sigma = 3 a): mostly long waves that reach the layer.
     * A one-node pluck also excites near-zone-boundary modes whose group
     * velocity is ~0; those linger in the interior whatever the layer does. */
    for (int r = 0; r < 64; r++)
        for (int c = 0; c < 64; c++) {
            double d2 = (c - 31.5) * (c - 31.5) + (r - 31.5) * (r - 31.5);
            s.uy[sim_idx(&s, c, r)] = 0.05 * exp(-d2 / 18.0);
        }
    sim_compute_forces(&s);
    double ke, pe; sim_energy(&s, &ke, &pe);
    double E0 = ke + pe;
    while (s.t < 6.0) sim_step(&s);
    sim_energy(&s, &ke, &pe);
    sim_free(&s);
    return (ke + pe) / E0;
}

static void test_absorb(void)
{
    printf("absorbing layer: energy left 6 s after a smooth central pulse\n");
    double off = absorb_residual(0, 0);
    double on  = absorb_residual(1, sqrt(400.0));
    CHECK(off > 0.99 && on < 0.01, "no layer %.3f, with layer %.4f", off, on);
}

/* 7a. crystal model: walls hold, cells never invert, thermostat is correct */
#ifndef KICKAMP
#define KICKAMP 25.0
#endif
static int any_inverted(Sim *s)
{
    int keep = s->p.orient, kf = s->p.bounded;
    s->p.orient = 1; s->p.bounded = 0;
    int ok = sim_valid(s);
    s->p.orient = keep; s->p.bounded = kf;
    return !ok;
}

static void test_crystal(void)
{
    printf("crystal model (bounded bonds, orientation barrier, thermostat)\n");

    /* energy conservation with bounded bonds at moderate amplitude */
    {
        Params p = base(32, BC_FREE, 1);
        p.bounded = 1;
        Sim s; sim_init(&s, &p);
        sim_kick_random(&s, 4.0, 4.0);
        double ke, pe; sim_energy(&s, &ke, &pe);
        double E0 = ke + pe, worst = 0;
        for (int i = 0; i < 40000; i++) {
            sim_step(&s);
            if (i % 100 == 0) { sim_energy(&s, &ke, &pe); double d = fabs(ke + pe - E0) / E0; if (d > worst) worst = d; }
        }
        CHECK(worst < 3e-4, "bounded, gamma = 0: max |dE|/E = %.1e over 20 s", worst);
        sim_free(&s);
    }

    /* violent kicks: without the barriers cells fold over; with them nothing
     * crosses a wall or inverts, and energy is still conserved */
    for (int guard = 0; guard < 2; guard++) {
        Params p = base(24, BC_PINNED, 1);
        p.bounded = guard; p.orient = guard;
        Sim s; sim_init(&s, &p);
        sim_kick_random(&s, KICKAMP, KICKAMP);
        double ke, pe; sim_energy(&s, &ke, &pe);
        double E0 = ke + pe, worst = 0;
        int inverted = 0, invalid = 0, maxdepth = 0, blocked = 0;
        for (int i = 0; i < 4000; i++) {
            sim_step(&s);
            if (s.sub_depth > maxdepth) maxdepth = s.sub_depth;
            blocked += s.rod_blocked;
            if (i % 10 == 0) {
                inverted += any_inverted(&s);
                if (guard && !sim_valid(&s)) invalid++;
                sim_energy(&s, &ke, &pe);
                double d = fabs(ke + pe - E0) / E0; if (d > worst) worst = d;
            }
        }
        if (!guard)
            CHECK(inverted > 0, "control (Hooke, no barrier, kick %.0f a/s): %d of 400 samples have folded cells",
                  KICKAMP, inverted);
        else
            CHECK(inverted == 0 && invalid == 0 && blocked == 0 && worst < 1e-2,
                  "bounded + barrier: 0 folds, 0 wall crossings; max substeps %d; reversals %ld; |dE|/E %.1e",
                  1 << maxdepth, s.reversals, worst);
        sim_free(&s);
    }

    /* pluck is capped inside the bond bound */
    {
        Params p = base(16, BC_FREE, 1);
        p.bounded = 1;
        Sim s; sim_init(&s, &p);
        int r = sim_displace(&s, sim_idx(&s, 8, 8), 0.0, -2.0);
        double d = fabs(s.uy[sim_idx(&s, 8, 8)]);
        CHECK(r == 0 && fabs(d - 0.8 * p.bound_delta) < 1e-12, "pluck of 2a capped to %.3f (0.8 delta)", d);
        sim_free(&s);
    }

    /* barrier alone refuses a pluck that would fold a cell */
    {
        Params p = base(16, BC_FREE, 1);
        p.orient = 1;
        Sim s; sim_init(&s, &p);
        int r = sim_displace(&s, sim_idx(&s, 8, 8), 0.0, -2.0);
        double u = 0; for (int i = 0; i < s.n; i++) u += fabs(s.uy[i]);
        CHECK(r == -1 && u == 0.0, "barrier: a 2a pluck (which folds cells) is refused");
        sim_free(&s);
    }

    /* thermostat: equipartition <KE> = kT per node (2 d.o.f.), and the
     * Lindemann ratio scales as sqrt(T) */
    double lind[2];
    for (int t = 0; t < 2; t++) {
        Params p = base(16, BC_PERIODIC, 1);
        p.thermo = 1; p.temp = t ? 2.0 : 0.5; p.therm_gamma = 4.0;
        Sim s; sim_init(&s, &p);
        double ksum = 0, l2 = 0, v2 = 0, v4 = 0; int ns = 0, nv = 0;
        for (int i = 0; i < 80000; i++) {
            sim_step(&s);
            if (i > 20000 && i % 20 == 0) {
                double ke, pe; sim_energy(&s, &ke, &pe);
                ksum += ke / s.n; double l = sim_lindemann(&s); l2 += l * l; ns++;
                for (int j = 0; j < s.n; j++) { double v = s.vx[j] * s.vx[j]; v2 += v; v4 += v * v; nv++; }
            }
        }
        double kT = ksum / ns, kurt = (v4 / nv) / ((v2 / nv) * (v2 / nv));
        lind[t] = sqrt(l2 / ns);
        CHECK(fabs(kT / p.temp - 1) < 0.05 && fabs(kurt - 3.0) < 0.1,
              "T = %.1f: <KE>/node = %.3f (expect %.3f), velocity kurtosis %.3f (Gaussian 3), Lindemann %.4f",
              p.temp, kT, p.temp, kurt, lind[t]);
        sim_free(&s);
    }
    CHECK(fabs(lind[1] / lind[0] - 2.0) < 0.1, "Lindemann ratio T=2 / T=0.5 = %.3f (expect 2 = sqrt 4)", lind[1] / lind[0]);
}

/* 7b. random kicks carry no net momentum or angular momentum (free) */
static void test_kick(void)
{
    printf("random kick: no net momentum / angular momentum (free BC)\n");
    Params p = base(16, BC_FREE, 1);
    Sim s; sim_init(&s, &p);
    sim_displace(&s, sim_idx(&s, 3, 5), 0.3, 0.2);
    sim_kick_random(&s, 5.0, 5.0);
    double px = 0, py = 0, Lz = 0, cx = 0, cy = 0, M = 0, vs = 0;
    for (int i = 0; i < s.n; i++) { M += s.m[i]; cx += s.m[i] * (i % 16 + s.ux[i]); cy += s.m[i] * (i / 16 + s.uy[i]); }
    cx /= M; cy /= M;
    for (int i = 0; i < s.n; i++) {
        double x = i % 16 + s.ux[i] - cx, y = i / 16 + s.uy[i] - cy;
        px += s.m[i] * s.vx[i]; py += s.m[i] * s.vy[i]; Lz += s.m[i] * (x * s.vy[i] - y * s.vx[i]);
        vs += fabs(s.vx[i]) + fabs(s.vy[i]);
    }
    CHECK(fabs(px) + fabs(py) < 1e-9 * vs && fabs(Lz) < 1e-9 * vs * 16, "P = (%.1e, %.1e), Lz = %.1e", px, py, Lz);
    sim_free(&s);
}

/* 8. throw-velocity fit */
static void test_throw(void)
{
    printf("release-velocity fit\n");
    MotionLog m; mlog_clear(&m);
    for (int i = 0; i < 40; i++) { double t = i / 60.0; mlog_push(&m, t, 3 + 2 * t, -1 - 5 * t); }
    double vx, vy; int n = mlog_velocity(&m, 39 / 60.0, 0.1, &vx, &vy);
    CHECK(fabs(vx - 2) < 1e-9 && fabs(vy + 5) < 1e-9 && n >= 2, "fit (%.6f, %.6f) from %d samples", vx, vy, n);
    mlog_clear(&m); mlog_push(&m, 0.0, 1, 1); mlog_push(&m, 0.01, 1, 1);
    n = mlog_velocity(&m, 0.5, 0.1, &vx, &vy);
    CHECK(vx == 0 && vy == 0, "stale samples give zero velocity (%d used)", n);
}

int main(int argc, char **argv)
{
    if (argc > 1 && argv[1][0] == 'g') {          /* ./test g : absorb gamma sweep */
        double off = absorb_residual(0, 0);
        printf("no layer: %.4f\n", off);
        for (double g = 5; g <= 160; g *= 2) printf("gamma %6.1f: %.4f\n", g, absorb_residual(1, g));
        return 0;
    }
    test_dynmat();
    test_energy();
    test_damping();
    test_dispersion();
    test_rods();
    test_mode_purity();
    test_fput();
    test_defect();
    test_absorb();
    test_crystal();
    test_kick();
    test_throw();
    printf(failures ? "\n%d FAILED\n" : "\nall passed\n", failures);
    return failures ? 1 : 0;
}
