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
    int     nspr, spr_cap;

    double prestrain_eff;  /* prestrain actually applied (0 unless clamped/periodic) */

    int    held;           /* dragged node, or -1 */
    double held_ux, held_uy, held_vx, held_vy;

    double    t;
    long long steps;
    double    pe_ref;      /* potential energy of the reference state */
    double    omega_max;   /* Gershgorin bound on the top frequency   */
    uint64_t  rng;
} Sim;

int    sim_init(Sim *s, const Params *p);
void   sim_free(Sim *s);

void   sim_set_masses(Sim *s);     /* after changing p.disorder / p.mass       */
void   sim_rebuild(Sim *s);        /* after bc, diagonal, prestrain, k changes */
void   sim_update_damping(Sim *s); /* after damping / absorb changes           */
void   sim_reset(Sim *s);          /* zero u and v                             */

void   sim_compute_forces(Sim *s);
void   sim_step(Sim *s);
void   sim_energy(const Sim *s, double *ke, double *pe);   /* pe relative to pe_ref */
void   sim_lattice(const Sim *s, Lattice *L);
double sim_mean_mass(const Sim *s);
double sim_min_mass(const Sim *s);
double sim_omega0(const Sim *s);   /* sqrt(k/m), a natural velocity/frequency scale */

void   sim_kick_random(Sim *s, double ampx, double ampy);   /* velocity kicks */
void   sim_displace(Sim *s, int i, double dx, double dy);
int    sim_nearest(const Sim *s, double x, double y, double maxd);
int    sim_defect_index(const Sim *s);                      /* -1 if none */
/* Replace the state with the lowest nonzero x-wavenumber mode of branch
 * b (0 = L, 1 = T) on the lowest y profile, peak strain `strain`. */
void   sim_excite_mode(Sim *s, int branch, double strain);

static inline int sim_idx(const Sim *s, int c, int r) { return r * s->cols + c; }

#endif
