#include "util.h"

#include <stdlib.h>
#include <string.h>

double rng_uniform(uint64_t *s)
{
    uint64_t z = (*s += 0x9E3779B97F4A7C15ull);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ull;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBull;
    z ^= z >> 31;
    return (double)(z >> 11) * (1.0 / 9007199254740992.0);
}

double rng_sym(uint64_t *s) { return 2.0 * rng_uniform(s) - 1.0; }

int hist_init(Hist *h, int nser, int cap, double dt)
{
    memset(h, 0, sizeof *h);
    h->nser = nser; h->cap = cap & ~1; h->dt0 = dt;
    h->t = calloc((size_t)h->cap, sizeof *h->t);
    h->v = calloc((size_t)h->cap * nser, sizeof *h->v);
    hist_clear(h, 0.0);
    return (h->t && h->v) ? 0 : -1;
}

void hist_free(Hist *h) { free(h->t); free(h->v); h->t = h->v = NULL; }

void hist_clear(Hist *h, double t_now)
{
    h->len = 0; h->dt = h->dt0; h->next_t = t_now; h->t0 = t_now;
}

void hist_push(Hist *h, double t, const double *vals)
{
    if (h->len == h->cap) {
        int k = 0;
        for (int i = 0; i < h->len; i += 2, k++) {
            h->t[k] = h->t[i];
            memcpy(&h->v[(size_t)k * h->nser], &h->v[(size_t)i * h->nser],
                   (size_t)h->nser * sizeof *h->v);
        }
        h->len = k;
        h->dt *= 2.0;
    }
    h->t[h->len] = t;
    memcpy(&h->v[(size_t)h->len * h->nser], vals, (size_t)h->nser * sizeof *vals);
    h->len++;
    h->next_t = t + h->dt;
}

void mlog_clear(MotionLog *m) { m->n = 0; m->head = 0; }

void mlog_push(MotionLog *m, double t, double x, double y)
{
    m->t[m->head] = t; m->x[m->head] = x; m->y[m->head] = y;
    m->head = (m->head + 1) % MLOG_N;
    if (m->n < MLOG_N) m->n++;
}

int mlog_velocity(const MotionLog *m, double t_now, double window, double *vx, double *vy)
{
    double st = 0, sx = 0, sy = 0, stt = 0, stx = 0, sty = 0;
    int k = 0;
    for (int i = 0; i < m->n; i++) {
        int j = (m->head - 1 - i + MLOG_N) % MLOG_N;
        double dt = m->t[j] - t_now;          /* <= 0 */
        if (dt < -window) break;
        st += dt; sx += m->x[j]; sy += m->y[j];
        stt += dt * dt; stx += dt * m->x[j]; sty += dt * m->y[j];
        k++;
    }
    *vx = *vy = 0.0;
    if (k < 2) return k;
    double den = k * stt - st * st;
    if (den <= 1e-18) return k;
    *vx = (k * stx - st * sx) / den;
    *vy = (k * sty - st * sy) / den;
    return k;
}
