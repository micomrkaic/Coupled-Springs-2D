/*
 * modes.h — normal-mode bases and the analytic dispersion relation
 *
 * 1-D bases along one lattice direction of N sites, one per boundary
 * condition, indexed by an integer wavenumber w:
 *
 *   free, pinned  DCT-II   phi_w(j) ~ cos(pi w (j+1/2) / N)   q = pi w / N,      w = 0..N-1
 *   clamped       DST-I    phi_w(j) ~ sin(pi w j / (N-1))     q = pi w / (N-1),  w = 1..N-2
 *   periodic      Fourier  cos and sin(2 pi w j / N)          q = 2 pi w / N,    w = 0..N/2
 *
 * All bases are orthonormal over the lattice sites.  For the free lattice
 * with cardinal springs only, DCT-II vectors are the exact eigenmodes; for
 * clamped, DST-I; for periodic, Fourier.  "Pinned" has no exact basis
 * and uses DCT-II as an approximation.
 *
 * Dynamical matrix of the infinite lattice, linearised about the
 * (possibly prestressed) reference state:
 *
 *   D(q) = (2/m) sum_{bonds r} K_r (1 - cos q.r),
 *   K_r  = k_r n n^T + (T_r/|r|)(I - n n^T),   T_r = k_r(|r| - L0_r)
 *
 * The eigenvalues of D are omega^2.  The eigenvector with the larger
 * overlap with q-hat is labelled quasi-longitudinal (L), the other quasi-
 * transverse (T).  Without diagonal springs or tension the eigenvectors
 * are x-hat and y-hat for every q, and omega_T = 0 along the axes.
 *
 * Rigid rods: a rod family with bond r imposes (e^{iq.r} - 1) n.u~ = 0,
 * i.e. n.u~ = 0 whenever q.r is not a multiple of 2 pi.  If the active
 * constraints leave one allowed direction t, that branch has
 * omega^2 = t^T D t and the other is frozen (omega = +inf); if they span
 * the plane, both are frozen.
 */
#ifndef MODES_H
#define MODES_H

#include "params.h"

int    line_w0(Boundary bc);                       /* lowest wavenumber    */
int    line_nw(Boundary bc, int N);                /* one past the highest */
double line_q(Boundary bc, int N, int w);
int    line_has_sin(Boundary bc, int N, int w);
/* s = 0: cosine (or the only) vector; s = 1: sine partner (periodic only) */
double line_phi(Boundary bc, int N, int w, int s, int j);

typedef struct {
    double k, kd;      /* cardinal and diagonal stiffness            */
    double tc, td;     /* tension / bond length, cardinal & diagonal */
    double m;          /* mass used for the analytic curves          */
    int    diag;
    int    rod_card, rod_diag;
} Lattice;

void dynmat(const Lattice *L, double qx, double qy, double D[3]);   /* xx, yy, xy */
/* w[0], e[0] = quasi-L; w[1], e[1] = quasi-T = e[0] rotated by +90 deg;
 * a branch frozen by rods has w = INFINITY */
void branches(const Lattice *L, double qx, double qy, double w[2], double e[2][2]);
double band_top(const Lattice *L);                 /* max omega over q */

#endif
