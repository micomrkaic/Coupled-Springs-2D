/*
 * util.h — small helpers with no SDL or FFTW dependency
 */
#ifndef UTIL_H
#define UTIL_H

#include <stdint.h>

/* splitmix64 */
double rng_uniform(uint64_t *state);              /* [0, 1)  */
double rng_sym(uint64_t *state);                  /* [-1, 1) */

/*
 * Decimating history: fixed capacity, sampled in simulation time.  When
 * full, every other sample is dropped and the sampling interval doubles,
 * so the whole run since the last clear stays on screen.
 */
typedef struct {
    int     nser, cap, len;
    double  dt0, dt, next_t, t0;
    double *t;          /* cap          */
    double *v;          /* cap * nser   */
} Hist;

int  hist_init(Hist *h, int nser, int cap, double dt);
void hist_free(Hist *h);
void hist_clear(Hist *h, double t_now);
void hist_push(Hist *h, double t, const double *vals);
static inline int hist_due(const Hist *h, double t) { return t >= h->next_t; }

/*
 * Pointer-motion log used to estimate the release ("throw") velocity of a
 * dragged mass by a least-squares fit of position on time over a short
 * window before release.
 */
#define MLOG_N 32
typedef struct { double t[MLOG_N], x[MLOG_N], y[MLOG_N]; int n, head; } MotionLog;

void mlog_clear(MotionLog *m);
void mlog_push(MotionLog *m, double t, double x, double y);
/* returns the number of samples used; vx = vy = 0 if fewer than 2 */
int  mlog_velocity(const MotionLog *m, double t_now, double window, double *vx, double *vy);

#endif
