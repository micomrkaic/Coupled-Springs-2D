/*
 * sim.h — spring-mass lattice: state, boundary conditions, integrator
 *
 * State is stored as displacements u from the reference lattice (site
 * (c, r) sits at (c, r) in lattice units), so periodic wrap-around costs
 * nothing: each spring carries its reference bond vector r and its current
 * extension is r + u_b - u_a.
 *
 * Integrator: Strang splitting
 *     damp(dt/2) . velocity-Verlet(dt) . damp(dt/2)
 * with the exact drag factor v <- v exp(-gamma_i dt / (2 m_i)).  The
 * conservative part is symplectic and time-reversible; the drag is solved
 * exactly, so gamma = 0 recovers plain Verlet.
 *
 * Rigid rods (per family: cardinal and/or diagonal) are holonomic
 * constraints |r + u_b - u_a| = |r|, enforced with RATTLE: SHAKE on the
 * positions after the drift, a velocity projection after the final kick.
 * RATTLE keeps the scheme symplectic.  Rods store no energy and do not
 * limit dt; their tension (the Lagrange multiplier) is kept for display.
 *
 * Crystal model (each switchable):
 *   bounded bonds  V(e) = 1/2 k e^2 / (1 - (e/D)^2),  e = r - L0: exactly
 *                  Hookean at small e, infinite at |e| = D, so a bond can
 *                  neither collapse (D < L0) nor stretch without bound.  The
 *                  wall rises like 1/distance; a logarithmic (FENE) wall is
 *                  too soft: at a few hundred k D^2 of energy its turning
 *                  point lies closer to the wall than double precision.
 *   orientation    every unit cell is split into its four corner triangles;
 *   barrier        each signed area A carries b(A) = kappa (A_on - A)^2 /
 *                  (A A_on) for A < A_on (zero above), which is C^1 at A_on
 *                  and infinite at A = 0, so no cell can fold over and
 *                  neighbours can never change their left/right order.
 *   thermostat     Langevin: the drag half-steps become exact Ornstein-
 *                  Uhlenbeck steps v <- c v + sqrt(kT/m (1 - c^2)) xi (OBABO).
 *   safeguard      a step that would cross a wall, or come so close that it
 *                  is not resolved (a bond closing more than half its gap to
 *                  the bound, a small triangle losing more than half its
 *                  area), or that SHAKE cannot close, is undone and redone as
 *                  2, 4, ... 1024 substeps.  If even that
 *                  fails, a drag is refused (rod_blocked); free motion is
 *                  time-reversed (v -> -v, energy-exact), counted in reversals.
 */
#ifndef SIM_H
#define SIM_H

#include "modes.h"
#include "params.h"

#include <stdint.h>

typedef struct {
    int    a, b;
    double rx, ry;      /* reference bond vector a -> b */
    double k, L0;
    int    wrap;        /* periodic wrap-around bond (not drawn) */
    int    rod;         /* rigid rod (constraint) instead of a spring */
    double f;           /* rod tension from the last step (N)          */
} Spring;

typedef struct {
    Params p;           /* live copy; runtime toggles edit this */
    int    cols, rows, n;

    double *ux, *uy, *vx, *vy, *ax, *ay;
    double *m, *inv_m;
    double *gamma;      /* total drag per node (uniform + absorbing layer) */
    double *decay;      /* exp(-gamma dt / 2m)                             */
    double *anchor_k;   /* on-site anchor stiffness (pinned BC)            */
    unsigned char *fixed;

    Spring *spr;
    int     nspr, spr_cap, nrods;
    double *oux, *ouy;                    /* positions at step start (SHAKE) */
    double *sx, *sy, *svx, *svy;          /* snapshot for undoing a blocked perturbation */
    int    *chain, *chain_off, nchains, max_chain;   /* rods grouped into straight chains */
    double *tri;                          /* tridiagonal scratch             */
    double *vfac;                         /* factored chain matrices (velocity solve) */
    int     rod_iters;                    /* SHAKE sweeps in the last step   */
    double  rod_resid;                    /* worst |d^2 - L^2| / L^2 after it */
    int     rod_blocked;                  /* last step was undone: constraints forbid the motion */
    int     sub_depth;                    /* last step used 2^sub_depth substeps */
    long    reversals;                    /* last-resort time reversals so far */
    double *sax, *say;                    /* accelerations at step start     */
    double *noise;                        /* thermostat noise amplitude per half step */
    int    *cell, ncells;                 /* unit cells: 4 corner indices each */
    double  lind_ema;                     /* smoothed Lindemann ratio (UI)   */
    double *exc;                          /* bond excursion r - L0 at the last force evaluation */
    double *tarea;                        /* triangle areas at the last force evaluation */
    double  wall_w2;                      /* (local wall frequency)^2 at the last force evaluation */

    double prestrain_eff;  /* prestrain actually applied (0 unless clamped/periodic) */

    int    held;           /* dragged node, or -1 */
    double held_ux, held_uy, held_vx, held_vy;   /* drag target (displacement, velocity) */
    int    drag_soft;      /* 0: held node follows the target exactly (kinematic);
                              1: a stiff, critically damped spring pulls it toward the
                                 target, with the pull capped at one lattice spacing, so
                                 walls and rods always win and nothing is refused */

    double    t;
    long long steps;
    double    pe_ref;      /* potential energy of the reference state */
    double    omega_max;   /* Gershgorin bound on the top frequency   */
    uint64_t  rng;
} Sim;

int    sim_init(Sim *s, const Params *p);
void   sim_free(Sim *s);

void   sim_set_masses(Sim *s);     /* after changing p.disorder / p.mass       */
/* after bc, diagonal, rods, prestrain, k changes; returns 1 if the state
 * had to be reset to the reference lattice to satisfy the rods */
int    sim_rebuild(Sim *s);
void   sim_update_damping(Sim *s); /* after damping / absorb changes           */
void   sim_reset(Sim *s);          /* zero u and v                             */

void   sim_compute_forces(Sim *s);
void   sim_step(Sim *s);
void   sim_energy(const Sim *s, double *ke, double *pe);   /* pe relative to pe_ref */
void   sim_lattice(const Sim *s, Lattice *L);
double sim_mean_mass(const Sim *s);
double sim_min_mass(const Sim *s);
double sim_omega0(const Sim *s);   /* sqrt(k/m), a natural velocity/frequency scale */
/* r.m.s. relative displacement of nearest neighbours, in units of a: the
 * 2-D Lindemann ratio (absolute displacements diverge in 2-D, Mermin-Wagner);
 * crystals melt near 0.1 */
double sim_lindemann(const Sim *s);
/* 1 if every bond is inside its bound and every cell keeps its
 * orientation (only the enabled checks are made) */
int    sim_valid(const Sim *s);

/* Perturbations.  With rods present each is projected onto the constraint
 * manifold; a displacement the rods cannot accommodate (e.g. clamped edges
 * with cardinal rods) is undone and returns -1. */
void   sim_kick_random(Sim *s, double ampx, double ampy);   /* velocity kicks */
int    sim_displace(Sim *s, int i, double dx, double dy);
/* project positions (nonlinear) and velocities onto the rod constraints;
 * returns the final relative residual (0 if no rods) */
double sim_project(Sim *s);
void   sim_project_velocities(Sim *s);
int    sim_nearest(const Sim *s, double x, double y, double maxd);
int    sim_defect_index(const Sim *s);                      /* -1 if none */
/* Replace the state with the lowest nonzero x-wavenumber mode of branch
 * b (0 = L, 1 = T) on the lowest y profile, peak strain `strain`. */
int    sim_excite_mode(Sim *s, int branch, double strain);

static inline int sim_idx(const Sim *s, int c, int r) { return r * s->cols + c; }

#endif
