#include "spectrum.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

#define DECADES_OCC  5.0
#define DECADES_DISP 4.0

static float *g_img = NULL;      /* scratch for the image getters */
static size_t g_img_cap = 0;

static float *img_buf(size_t n)
{
    if (n > g_img_cap) {
        float *p = realloc(g_img, n * sizeof *p);
        if (!p) return NULL;
        g_img = p; g_img_cap = n;
    }
    return g_img;
}

/* ───────────────────────────────────────────────────────── setup ── */

static double q_axis(Boundary bc, int N, int k, int full_signed)
{
    switch (bc) {
    case BC_PERIODIC: {
        int kk = (full_signed && k > N / 2) ? k - N : k;
        return 2.0 * PI * kk / N;
    }
    case BC_CLAMPED: return PI * (k + 1) / (N - 1);
    default:         return PI * k / N;
    }
}

int spectrum_init(Spectrum *sp, const Sim *s)
{
    memset(sp, 0, sizeof *sp);
    const Params *p = &s->p;
    sp->bc = p->bc; sp->rows = s->rows; sp->cols = s->cols;
    sp->periodic = p->bc == BC_PERIODIC;
    sp->off = p->bc == BC_CLAMPED ? 1 : 0;
    sp->ny = s->rows - 2 * sp->off;
    sp->nx = s->cols - 2 * sp->off;
    sp->nqy = sp->ny;
    sp->nqx = sp->periodic ? sp->nx / 2 + 1 : sp->nx;
    sp->occ_alpha = 0.15;

    size_t nin = (size_t)sp->ny * sp->nx, nq = (size_t)sp->nqy * sp->nqx;
    sp->rin = fftw_alloc_real(nin);
    if (sp->periodic) {
        sp->cout = fftw_alloc_complex(nq);
        sp->plan = fftw_plan_dft_r2c_2d(sp->ny, sp->nx, sp->rin, sp->cout, FFTW_ESTIMATE);
    } else {
        fftw_r2r_kind kind = p->bc == BC_CLAMPED ? FFTW_RODFT00 : FFTW_REDFT10;
        sp->rout = fftw_alloc_real(nin);
        sp->plan = fftw_plan_r2r_2d(sp->ny, sp->nx, sp->rin, sp->rout, kind, kind, FFTW_ESTIMATE);
    }
    if (!sp->rin || !sp->plan) return -1;

    double **arrs[] = { &sp->Xr, &sp->Xi, &sp->Yr, &sp->Yi, &sp->wL, &sp->wT,
                        &sp->eLx, &sp->eLy, &sp->occL, &sp->occT };
    for (size_t i = 0; i < sizeof arrs / sizeof *arrs; i++)
        if (!(*arrs[i] = calloc(nq, sizeof(double)))) return -1;
    sp->qx = calloc((size_t)sp->nqx, sizeof(double));
    sp->qy = calloc((size_t)sp->nqy, sizeof(double));
    if (!sp->qx || !sp->qy) return -1;
    for (int ix = 0; ix < sp->nqx; ix++) sp->qx[ix] = q_axis(p->bc, s->cols, ix, 0);
    for (int iy = 0; iy < sp->nqy; iy++) sp->qy[iy] = q_axis(p->bc, s->rows, iy, 1);

    sim_lattice(s, &sp->lat);
    for (int iy = 0; iy < sp->nqy; iy++)
        for (int ix = 0; ix < sp->nqx; ix++) {
            double w[2], e[2][2];
            branches(&sp->lat, sp->qx[ix], sp->qy[iy], w, e);
            size_t k = (size_t)iy * sp->nqx + ix;
            sp->wL[k] = w[0]; sp->wT[k] = w[1];
            sp->eLx[k] = e[0][0]; sp->eLy[k] = e[0][1];
        }

    /* omega-q sampling: Nyquist ~1.25x the top of the band (or of a light
     * defect's local mode, which sits above the band) */
    double top = band_top(&sp->lat);
    double kxx = 2 * sp->lat.k + 2 * sp->lat.tc + (sp->lat.diag ? 2 * (sp->lat.kd + sp->lat.td) : 0);
    double wloc = sqrt(kxx / sim_min_mass(s)) * 1.05;
    if (wloc > top) top = wloc;
    double dt_target = PI / (1.25 * top);
    sp->stride = (int)floor(dt_target / p->dt);
    if (sp->stride < 1) sp->stride = 1;
    sp->dt_s = sp->stride * p->dt;
    sp->omega_nyq = PI / sp->dt_s;

    sp->nt = p->nt;
    sp->nlq = sp->nqx;
    size_t nl = (size_t)sp->nt * sp->nlq, nh = (size_t)(sp->nt / 2 + 1) * sp->nlq;
    for (int b = 0; b < 2; b++) {
        sp->lbuf[b] = fftw_alloc_complex(nl);
        sp->dpow[b] = calloc(nh, sizeof(double));
        if (!sp->lbuf[b] || !sp->dpow[b]) return -1;
        memset(sp->lbuf[b], 0, nl * sizeof(fftw_complex));
    }
    sp->tbuf = fftw_alloc_complex(nl);
    sp->hann = calloc((size_t)sp->nt, sizeof(double));
    if (!sp->tbuf || !sp->hann) return -1;
    int n = sp->nt;
    sp->tplan = fftw_plan_many_dft(1, &n, sp->nlq, sp->tbuf, NULL, 1, n,
                                   sp->tbuf, NULL, 1, n, FFTW_FORWARD, FFTW_ESTIMATE);
    if (!sp->tplan) return -1;

    /* mode-energy bases */
    sp->wx0 = 1;
    int avail = line_nw(p->bc, s->cols) - sp->wx0;
    sp->K = avail < SPEC_KMODES ? avail : SPEC_KMODES;
    sp->phix = calloc((size_t)sp->K * 2 * s->cols, sizeof(double));
    sp->psiy = calloc((size_t)s->rows, sizeof(double));
    if (!sp->phix || !sp->psiy) return -1;
    int wy = line_w0(p->bc);
    double qy0 = line_q(p->bc, s->rows, wy);
    for (int r = 0; r < s->rows; r++) sp->psiy[r] = line_phi(p->bc, s->rows, wy, 0, r);
    for (int k = 0; k < sp->K; k++) {
        int w = sp->wx0 + k;
        for (int si = 0; si < 2; si++)
            for (int c = 0; c < s->cols; c++)
                sp->phix[((size_t)k * 2 + si) * s->cols + c] =
                    (si && !line_has_sin(p->bc, s->cols, w)) ? 0.0 : line_phi(p->bc, s->cols, w, si, c);
        double wv[2], e[2][2];
        branches(&sp->lat, line_q(p->bc, s->cols, w), qy0, wv, e);
        sp->wmode[0][k] = wv[0]; sp->wmode[1][k] = wv[1];
    }
    if (hist_init(&sp->modes, 2 * SPEC_KMODES, 1024, 0.02) ||
        hist_init(&sp->energy, 3, 1024, 0.02)) return -1;

    spectrum_clear_history(sp, s);
    return 0;
}

void spectrum_free(Spectrum *sp)
{
    if (sp->plan)  fftw_destroy_plan(sp->plan);
    if (sp->tplan) fftw_destroy_plan(sp->tplan);
    fftw_free(sp->rin); fftw_free(sp->rout); fftw_free(sp->cout); fftw_free(sp->tbuf);
    for (int b = 0; b < 2; b++) { fftw_free(sp->lbuf[b]); free(sp->dpow[b]); }
    free(sp->Xr); free(sp->Xi); free(sp->Yr); free(sp->Yi);
    free(sp->wL); free(sp->wT); free(sp->eLx); free(sp->eLy);
    free(sp->occL); free(sp->occT); free(sp->qx); free(sp->qy);
    free(sp->hann); free(sp->phix); free(sp->psiy);
    hist_free(&sp->modes); hist_free(&sp->energy);
    memset(sp, 0, sizeof *sp);
}

int spectrum_rebuild(Spectrum *sp, const Sim *s)
{
    spectrum_free(sp);
    return spectrum_init(sp, s);
}

void spectrum_clear_history(Spectrum *sp, const Sim *s)
{
    sp->nsamp = 0; sp->head = 0; sp->dirty = 1;
    size_t nq = (size_t)sp->nqy * sp->nqx;
    memset(sp->occL, 0, nq * sizeof(double));
    memset(sp->occT, 0, nq * sizeof(double));
    hist_clear(&sp->modes, s->t);
    hist_clear(&sp->energy, s->t);
}

/* ───────────────────────────────────────────────────────── transforms ── */

static void transform_field(Spectrum *sp, const double *f, int cols, double *Re, double *Im)
{
    for (int iy = 0; iy < sp->ny; iy++)
        for (int ix = 0; ix < sp->nx; ix++)
            sp->rin[iy * sp->nx + ix] = f[(iy + sp->off) * cols + ix + sp->off];
    fftw_execute(sp->plan);
    size_t nq = (size_t)sp->nqy * sp->nqx;
    if (sp->periodic) {
        for (size_t k = 0; k < nq; k++) { Re[k] = sp->cout[k][0]; Im[k] = sp->cout[k][1]; }
    } else {
        memcpy(Re, sp->rout, nq * sizeof(double));
        memset(Im, 0, nq * sizeof(double));
    }
}

static void transform_velocity(Spectrum *sp, const Sim *s)
{
    transform_field(sp, s->vx, s->cols, sp->Xr, sp->Xi);
    transform_field(sp, s->vy, s->cols, sp->Yr, sp->Yi);
}

/* project (X, Y) at index k onto L and T: complex amplitudes */
static inline void project(const Spectrum *sp, size_t k,
                           double *Lr, double *Li, double *Tr, double *Ti)
{
    double ex = sp->eLx[k], ey = sp->eLy[k];
    *Lr =  ex * sp->Xr[k] + ey * sp->Yr[k];
    *Li =  ex * sp->Xi[k] + ey * sp->Yi[k];
    *Tr = -ey * sp->Xr[k] + ex * sp->Yr[k];
    *Ti = -ey * sp->Xi[k] + ex * sp->Yi[k];
}

/* ───────────────────────────────────────────────────────── sampling ── */

static void sample_modes(Spectrum *sp, const Sim *s)
{
    const int C = s->cols, R = s->rows, K = sp->K;
    double Ux[1024], Uy[1024], Vx[1024], Vy[1024];
    for (int c = 0; c < C; c++) {
        double ux = 0, uy = 0, vx = 0, vy = 0;
        for (int r = 0; r < R; r++) {
            int i = r * C + c; double w = sp->psiy[r];
            ux += w * s->ux[i]; uy += w * s->uy[i];
            vx += w * s->vx[i]; vy += w * s->vy[i];
        }
        Ux[c] = ux; Uy[c] = uy; Vx[c] = vx; Vy[c] = vy;
    }
    double m = sp->lat.m, out[2 * SPEC_KMODES] = { 0 };
    int wy = line_w0(sp->bc);
    double qy0 = line_q(sp->bc, R, wy);
    for (int k = 0; k < K; k++) {
        double wv[2], e[2][2];
        branches(&sp->lat, line_q(sp->bc, C, sp->wx0 + k), qy0, wv, e);
        double EL = 0, ET = 0;
        for (int si = 0; si < 2; si++) {
            const double *ph = &sp->phix[((size_t)k * 2 + si) * C];
            double qx = 0, qy = 0, px = 0, py = 0;
            for (int c = 0; c < C; c++) {
                qx += ph[c] * Ux[c]; qy += ph[c] * Uy[c];
                px += ph[c] * Vx[c]; py += ph[c] * Vy[c];
            }
            double QL = e[0][0] * qx + e[0][1] * qy, PL = e[0][0] * px + e[0][1] * py;
            double QT = e[1][0] * qx + e[1][1] * qy, PT = e[1][0] * px + e[1][1] * py;
            if (isfinite(wv[0])) EL += 0.5 * m * (PL * PL + wv[0] * wv[0] * QL * QL);
            if (isfinite(wv[1])) ET += 0.5 * m * (PT * PT + wv[1] * wv[1] * QT * QT);
        }
        out[k] = EL; out[SPEC_KMODES + k] = ET;
    }
    hist_push(&sp->modes, s->t, out);
}

void spectrum_after_step(Spectrum *sp, const Sim *s)
{
    if (s->steps % sp->stride == 0) {
        transform_velocity(sp, s);
        size_t base = (size_t)sp->head * sp->nlq;
        for (int ix = 0; ix < sp->nlq; ix++) {
            double Lr, Li, Tr, Ti;
            project(sp, (size_t)ix, &Lr, &Li, &Tr, &Ti);     /* row iy = 0 */
            sp->lbuf[0][base + ix][0] = Lr; sp->lbuf[0][base + ix][1] = Li;
            sp->lbuf[1][base + ix][0] = Tr; sp->lbuf[1][base + ix][1] = Ti;
        }
        sp->head = (sp->head + 1) % sp->nt;
        if (sp->nsamp < sp->nt) sp->nsamp++;
        sp->dirty = 1;
    }
    if (hist_due(&sp->modes, s->t) && sp->K > 0 && s->cols <= 1024) sample_modes(sp, s);
    if (hist_due(&sp->energy, s->t)) {
        double ke, pe, v[3];
        sim_energy(s, &ke, &pe);
        v[0] = ke; v[1] = pe; v[2] = ke + pe;
        hist_push(&sp->energy, s->t, v);
    }
}

void spectrum_frame(Spectrum *sp, const Sim *s)
{
    transform_velocity(sp, s);
    size_t nq = (size_t)sp->nqy * sp->nqx;
    double a = sp->occ_alpha;
    for (size_t k = 0; k < nq; k++) {
        double Lr, Li, Tr, Ti;
        project(sp, k, &Lr, &Li, &Tr, &Ti);
        sp->occL[k] = (1 - a) * sp->occL[k] + a * (Lr * Lr + Li * Li);
        sp->occT[k] = (1 - a) * sp->occT[k] + a * (Tr * Tr + Ti * Ti);
    }
}

/* ───────────────────────────────────────────────────────── images ── */

static inline double occ_val(const Spectrum *sp, Pol pol, size_t k)
{
    return pol == POL_L ? sp->occL[k] : pol == POL_T ? sp->occT[k] : sp->occL[k] + sp->occT[k];
}

static inline float log_norm(double v, double mx, double decades)
{
    if (v <= 0 || mx <= 0) return 0.0f;
    double t = log10(v / mx) / decades + 1.0;
    return (float)(t < 0 ? 0 : t > 1 ? 1 : t);
}

const float *spectrum_occ_image(Spectrum *sp, Pol pol, int *w, int *h)
{
    int W, H;
    if (sp->periodic) { W = sp->nx; H = sp->ny; }
    else              { W = sp->nqx; H = sp->nqy; }
    float *img = img_buf((size_t)W * H);
    if (!img) return NULL;

    double mx = 0;
    size_t nq = (size_t)sp->nqy * sp->nqx;
    for (size_t k = 0; k < nq; k++) { double v = occ_val(sp, pol, k); if (v > mx) mx = v; }

    for (int dr = 0; dr < H; dr++)
        for (int dc = 0; dc < W; dc++) {
            size_t k;
            if (sp->periodic) {
                int ky = H / 2 - dr, kx = dc - W / 2;     /* qy up, qx right, origin centred */
                if (kx < 0) { kx = -kx; ky = -ky; }       /* |F(-q)| = |F(q)| */
                int iy = ((ky % H) + H) % H;
                k = (size_t)iy * sp->nqx + kx;
            } else {
                k = (size_t)(H - 1 - dr) * sp->nqx + dc;  /* origin bottom-left */
            }
            img[(size_t)dr * W + dc] = log_norm(occ_val(sp, pol, k), mx, DECADES_OCC);
        }
    *w = W; *h = H;
    return img;
}

static void compute_dispersion(Spectrum *sp)
{
    const int nt = sp->nt, nl = sp->nlq, nh = nt / 2 + 1;
    int L = sp->nsamp < nt ? sp->nsamp : nt;
    for (int b = 0; b < 2; b++) {
        if (L < 4) { memset(sp->dpow[b], 0, (size_t)nh * nl * sizeof(double)); continue; }
        for (int t = 0; t < L; t++) sp->hann[t] = 0.5 * (1.0 - cos(2.0 * PI * t / (L - 1)));
        for (int ix = 0; ix < nl; ix++) {
            fftw_complex *dst = &sp->tbuf[(size_t)ix * nt];
            for (int t = 0; t < nt; t++) {
                if (t < L) {
                    int slot = sp->nsamp < nt ? t : (sp->head + t) % nt;
                    dst[t][0] = sp->hann[t] * sp->lbuf[b][(size_t)slot * nl + ix][0];
                    dst[t][1] = sp->hann[t] * sp->lbuf[b][(size_t)slot * nl + ix][1];
                } else {
                    dst[t][0] = dst[t][1] = 0.0;
                }
            }
        }
        fftw_execute(sp->tplan);
        for (int ix = 0; ix < nl; ix++) {
            fftw_complex *A = &sp->tbuf[(size_t)ix * nt];
            for (int iw = 0; iw < nh; iw++) {
                double p = A[iw][0] * A[iw][0] + A[iw][1] * A[iw][1];
                if (iw > 0 && iw < nt / 2)
                    p += A[nt - iw][0] * A[nt - iw][0] + A[nt - iw][1] * A[nt - iw][1];
                sp->dpow[b][(size_t)iw * nl + ix] = p;
            }
        }
    }
    sp->dirty = 0;
}

const float *spectrum_disp_image(Spectrum *sp, Pol pol, int *w, int *h)
{
    if (sp->dirty) compute_dispersion(sp);
    const int W = sp->nlq, H = sp->nt / 2 + 1;
    float *img = img_buf((size_t)W * H);
    if (!img) return NULL;

    double gmax = 0;
    for (size_t k = 0; k < (size_t)W * H; k++) {
        double v = pol == POL_L ? sp->dpow[0][k] : pol == POL_T ? sp->dpow[1][k]
                                                 : sp->dpow[0][k] + sp->dpow[1][k];
        if (v > gmax) gmax = v;
    }
    /* normalise per q column (with a floor) so every column shows where in
     * omega its power sits */
    for (int ix = 0; ix < W; ix++) {
        double cmax = 0;
        for (int iw = 0; iw < H; iw++) {
            size_t k = (size_t)iw * W + ix;
            double v = pol == POL_L ? sp->dpow[0][k] : pol == POL_T ? sp->dpow[1][k]
                                                     : sp->dpow[0][k] + sp->dpow[1][k];
            if (v > cmax) cmax = v;
        }
        if (cmax < 1e-4 * gmax) cmax = 1e-4 * gmax;
        for (int iw = 0; iw < H; iw++) {
            size_t k = (size_t)iw * W + ix;
            double v = pol == POL_L ? sp->dpow[0][k] : pol == POL_T ? sp->dpow[1][k]
                                                     : sp->dpow[0][k] + sp->dpow[1][k];
            img[(size_t)(H - 1 - iw) * W + ix] = log_norm(v, cmax, DECADES_DISP);
        }
    }
    *w = W; *h = H;
    return img;
}
