/*
 * spectrum.h — spectral diagnostics of the lattice
 *
 * Spatial transform: matched to the boundary condition
 *     free / pinned   2-D DCT-II  (FFTW REDFT10) over all sites
 *     clamped         2-D DST-I   (FFTW RODFT00) over interior sites
 *     periodic        2-D real-to-complex FFT
 * and applied to both velocity components.  At each q the pair
 * (vx~, vy~) is projected onto the analytic branch eigenvectors e_L(q),
 * e_T(q) (see modes.h).  This split is exact for the periodic lattice and
 * approximate for the cosine/sine bases when diagonal springs are on.
 *
 * Views
 *   occupation   |v_b(q)|^2, exponentially averaged: which modes hold energy.
 *                It contains no frequency information.
 *   dispersion   |v_b(qx, omega)|^2 on the lowest-qy line: a space-time
 *                Fourier transform of the last nt samples (Hann window),
 *                with the analytic omega_L(q), omega_T(q) overlaid.
 *   modes        energy in the lowest K x-wavenumbers on the lowest-qy line,
 *                per branch, vs time (FPUT-style plot).
 *   energy       KE, PE, total vs time.
 */
#ifndef SPECTRUM_H
#define SPECTRUM_H

#include "sim.h"
#include "util.h"

#include <fftw3.h>

typedef enum { VIEW_OCC, VIEW_DISP, VIEW_MODES, VIEW_ENERGY, VIEW_COUNT } View;
typedef enum { POL_L, POL_T, POL_BOTH, POL_COUNT } Pol;

#define SPEC_KMODES 8

typedef struct {
    Boundary bc;
    int      rows, cols;
    int      ny, nx, off;      /* transform domain and its site offset       */
    int      nqy, nqx;         /* transform output dims                      */
    int      periodic;

    double       *rin, *rout;
    fftw_complex *cout;
    fftw_plan     plan;

    double *Xr, *Xi, *Yr, *Yi;   /* transformed vx, vy             */
    double *qx, *qy;             /* q for each output index         */
    double *wL, *wT, *eLx, *eLy; /* analytic branches per (iy, ix)  */
    double *occL, *occT;         /* EMA of |v_b(q)|^2               */
    double  occ_alpha;

    /* omega-q */
    int           nt, nlq, stride, nsamp, head, dirty;
    double        dt_s, omega_nyq;
    fftw_complex *lbuf[2];       /* nt * nlq ring                    */
    fftw_complex *tbuf;          /* nlq * nt                         */
    fftw_plan     tplan;
    double       *hann;
    double       *dpow[2];       /* (nt/2+1) * nlq                   */

    /* mode energies and total energy */
    int     K;                   /* number of x-wavenumbers shown   */
    int     wx0;                 /* first wavenumber shown          */
    double *phix;                /* K*2 * cols basis vectors        */
    double *psiy;                /* rows                            */
    double  wmode[2][SPEC_KMODES];
    Hist    modes;               /* 2K series: L then T             */
    Hist    energy;              /* KE, PE, E                       */

    Lattice lat;
} Spectrum;

int  spectrum_init(Spectrum *sp, const Sim *s);
void spectrum_free(Spectrum *sp);
int  spectrum_rebuild(Spectrum *sp, const Sim *s);  /* after any structural change */
void spectrum_clear_history(Spectrum *sp, const Sim *s);

void spectrum_after_step(Spectrum *sp, const Sim *s);  /* call after each sim_step */
void spectrum_frame(Spectrum *sp, const Sim *s);       /* call once per rendered frame */

/* Images in display orientation, values in [0,1] (log scale).
 * Pointers stay valid until the next call. */
const float *spectrum_occ_image(Spectrum *sp, Pol pol, int *w, int *h);
const float *spectrum_disp_image(Spectrum *sp, Pol pol, int *w, int *h);

#endif
