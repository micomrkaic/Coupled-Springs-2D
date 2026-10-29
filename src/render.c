#include "render.h"
#include "font.h"

#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* ───────────────────────────────────────────────────────── setup ── */

int render_init(Renderer *R, SDL_Window *win, SDL_Renderer *ren)
{
    memset(R, 0, sizeof *R);
    R->win = win; R->ren = ren;
    R->vcap = 1 << 16; R->icap = 1 << 17;
    R->vtx = malloc((size_t)R->vcap * sizeof *R->vtx);
    R->idx = malloc((size_t)R->icap * sizeof *R->idx);
    return (R->vtx && R->idx) ? 0 : -1;
}

void render_free(Renderer *R)
{
    if (R->tex) SDL_DestroyTexture(R->tex);
    free(R->vtx); free(R->idx);
    memset(R, 0, sizeof *R);
}

void render_layout(Renderer *R, const Sim *s)
{
    SDL_GetRendererOutputSize(R->ren, &R->out_w, &R->out_h);
    SDL_GetWindowSize(R->win, &R->win_w, &R->win_h);
    R->dpi = R->win_w > 0 ? (float)R->out_w / R->win_w : 1.0f;
    R->fs  = (int)lroundf(R->dpi); if (R->fs < 1) R->fs = 1;

    int pad = 8 * R->fs, lh = font_line_height(R->fs);
    int pw = (int)(330 * R->dpi);
    if (pw > R->out_w * 45 / 100) pw = R->out_w * 45 / 100;
    R->sim = (SDL_Rect){ 0, 0, R->out_w - pw, R->out_h };

    int side = pw - 2 * pad;
    if (side > R->out_h * 45 / 100) side = R->out_h * 45 / 100;
    R->spec = (SDL_Rect){ R->sim.w + pad, pad + lh + 2 * R->fs, side, side };
    int ty = R->spec.y + side + 2 * lh;
    R->text = (SDL_Rect){ R->sim.w + pad, ty, pw - 2 * pad, R->out_h - ty - pad };

    double sx = (double)R->sim.w / (s->cols + 1), sy = (double)R->sim.h / (s->rows + 1);
    R->scale = sx < sy ? sx : sy;
    R->ox = R->sim.x + 0.5 * (R->sim.w - (s->cols - 1) * R->scale);
    R->oy = R->sim.y + 0.5 * (R->sim.h - (s->rows - 1) * R->scale);
}

int render_in_sim(const Renderer *R, int wx, int wy)
{
    int x = (int)(wx * R->dpi), y = (int)(wy * R->dpi);
    return x >= R->sim.x && x < R->sim.x + R->sim.w && y >= R->sim.y && y < R->sim.y + R->sim.h;
}

void render_to_lattice(const Renderer *R, int wx, int wy, double *lx, double *ly)
{
    *lx = (wx * R->dpi - R->ox) / R->scale;
    *ly = (wy * R->dpi - R->oy) / R->scale;
}

/* ───────────────────────────────────────────────────────── geometry batch ── */

static void geo_reset(Renderer *R) { R->nv = R->ni = 0; }

static int geo_reserve(Renderer *R, int nv, int ni)
{
    if (R->nv + nv > R->vcap) {
        int cap = R->vcap; while (cap < R->nv + nv) cap *= 2;
        SDL_Vertex *p = realloc(R->vtx, (size_t)cap * sizeof *p);
        if (!p) return -1;
        R->vtx = p; R->vcap = cap;
    }
    if (R->ni + ni > R->icap) {
        int cap = R->icap; while (cap < R->ni + ni) cap *= 2;
        int *p = realloc(R->idx, (size_t)cap * sizeof *p);
        if (!p) return -1;
        R->idx = p; R->icap = cap;
    }
    return 0;
}

static inline int geo_v(Renderer *R, float x, float y, SDL_Color c)
{
    R->vtx[R->nv] = (SDL_Vertex){ { x, y }, c, { 0, 0 } };
    return R->nv++;
}

static inline void geo_t(Renderer *R, int a, int b, int c)
{
    R->idx[R->ni++] = a; R->idx[R->ni++] = b; R->idx[R->ni++] = c;
}

static void geo_line(Renderer *R, float x0, float y0, float x1, float y1, float w, SDL_Color c)
{
    if (geo_reserve(R, 4, 6)) return;
    float dx = x1 - x0, dy = y1 - y0, len = sqrtf(dx * dx + dy * dy);
    if (len < 1e-3f) return;
    float nx = -dy / len * 0.5f * w, ny = dx / len * 0.5f * w;
    int a = geo_v(R, x0 + nx, y0 + ny, c), b = geo_v(R, x1 + nx, y1 + ny, c);
    int d = geo_v(R, x1 - nx, y1 - ny, c), e = geo_v(R, x0 - nx, y0 - ny, c);
    geo_t(R, a, b, d); geo_t(R, a, d, e);
}

static void geo_disc(Renderer *R, float cx, float cy, float r, SDL_Color c)
{
    enum { SEG = 8 };
    if (geo_reserve(R, SEG + 1, 3 * SEG)) return;
    int ctr = geo_v(R, cx, cy, c), first = R->nv;
    for (int k = 0; k < SEG; k++) {
        float th = (float)(2 * PI * (k + 0.5) / SEG);
        geo_v(R, cx + r * cosf(th), cy + r * sinf(th), c);
    }
    for (int k = 0; k < SEG; k++) geo_t(R, ctr, first + k, first + (k + 1) % SEG);
}

static void geo_rect(Renderer *R, float x, float y, float w, float h, SDL_Color c)
{
    if (geo_reserve(R, 4, 6)) return;
    int a = geo_v(R, x, y, c), b = geo_v(R, x + w, y, c);
    int d = geo_v(R, x + w, y + h, c), e = geo_v(R, x, y + h, c);
    geo_t(R, a, b, d); geo_t(R, a, d, e);
}

static void geo_flush(Renderer *R)
{
    if (R->ni) SDL_RenderGeometry(R->ren, NULL, R->vtx, R->nv, R->idx, R->ni);
    geo_reset(R);
}

/* ───────────────────────────────────────────────────────── colours ── */

static SDL_Color strain_colour(double strain)
{
    double t = fabs(strain) * 8.0; if (t > 1) t = 1;
    Uint8 v = (Uint8)(255 - 200 * t), a = (Uint8)(80 + 100 * t);
    return strain > 0 ? (SDL_Color){ 255, v, v, a } : (SDL_Color){ v, v, 255, a };
}

static SDL_Color speed_colour(double spd, double vref)
{
    double t = spd / vref; if (t > 1) t = 1;
    return (SDL_Color){ (Uint8)(50 + 205 * t), (Uint8)(130 - 80 * t), (Uint8)(220 - 200 * t), 255 };
}

static Uint32 heat(float t)
{
    Uint8 r, g, b;
    if (t < 0.5f) { float s = t * 2; r = (Uint8)(120 * s); g = 0; b = (Uint8)(180 * s); }
    else { float s = (t - 0.5f) * 2; r = (Uint8)(120 + 135 * s); g = (Uint8)(220 * s); b = (Uint8)(180 + 75 * s); }
    return 0xFF000000u | (Uint32)r << 16 | (Uint32)g << 8 | b;
}

static const SDL_Color MODE_COL[SPEC_KMODES] = {
    { 255,  90,  90, 255 }, { 255, 170,  60, 255 }, { 240, 230,  80, 255 }, { 120, 220,  90, 255 },
    {  70, 210, 210, 255 }, {  90, 140, 255, 255 }, { 180, 110, 255, 255 }, { 240, 120, 200, 255 },
};
static const SDL_Color COL_L = { 80, 220, 255, 255 }, COL_T = { 255, 210, 80, 255 };

static void set_col(SDL_Renderer *r, SDL_Color c) { SDL_SetRenderDrawColor(r, c.r, c.g, c.b, c.a); }

/* ───────────────────────────────────────────────────────── lattice ── */

static void render_lattice(Renderer *R, const Sim *s)
{
    SDL_Renderer *ren = R->ren;
    const double sc = R->scale;
    SDL_RenderSetClipRect(ren, &R->sim);
    geo_reset(R);

#define PX(i) ((float)(R->ox + ((i) % s->cols + s->ux[i]) * sc))
#define PY(i) ((float)(R->oy + ((i) / s->cols + s->uy[i]) * sc))

    /* absorbing layer */
    if (s->p.absorb) {
        int W = s->p.absorb_width;
        float x0 = (float)(R->ox - 0.5 * sc), y0 = (float)(R->oy - 0.5 * sc);
        float w = (float)(s->cols * sc), h = (float)(s->rows * sc);
        for (int k = 0; k < W; k++) {      /* graded shading, darker outward */
            float f = (float)(W - k) / W, o = (float)(k * sc);
            SDL_Color c = { 150, 90, 40, (Uint8)(15 + 75 * f * f) };
            geo_rect(R, x0 + o, y0 + o, w - 2 * o, (float)sc, c);
            geo_rect(R, x0 + o, y0 + h - o - (float)sc, w - 2 * o, (float)sc, c);
            geo_rect(R, x0 + o, y0 + o + (float)sc, (float)sc, h - 2 * o - 2 * (float)sc, c);
            geo_rect(R, x0 + w - o - (float)sc, y0 + o + (float)sc, (float)sc, h - 2 * o - 2 * (float)sc, c);
        }
        geo_flush(R);
    }

    /* springs */
    float lw = (float)(sc * 0.08); if (lw < R->dpi) lw = R->dpi;
    for (int k = 0; k < s->nspr; k++) {
        const Spring *sp = &s->spr[k];
        if (sp->wrap) continue;
        int a = sp->a, b = sp->b;
        double dx = sp->rx + s->ux[b] - s->ux[a], dy = sp->ry + s->uy[b] - s->uy[a];
        double strain = (sqrt(dx * dx + dy * dy) - sp->L0) / sp->L0;
        geo_line(R, PX(a), PY(a), PX(b), PY(b), lw, strain_colour(strain));
    }
    geo_flush(R);

    /* anchors (pinned) and fixed nodes (clamped) */
    if (s->p.bc == BC_PINNED || s->p.bc == BC_CLAMPED) {
        for (int i = 0; i < s->n; i++) {
            float rx = (float)(R->ox + (i % s->cols) * sc), ry = (float)(R->oy + (i / s->cols) * sc);
            if (s->anchor_k[i] > 0) {
                geo_line(R, rx, ry, PX(i), PY(i), lw, (SDL_Color){ 80, 80, 200, 160 });
                geo_rect(R, rx - 1.5f * R->dpi, ry - 1.5f * R->dpi, 3 * R->dpi, 3 * R->dpi, (SDL_Color){ 60, 60, 160, 220 });
            }
            if (s->fixed[i]) {
                float h = (float)(sc * 0.22); if (h < 2 * R->dpi) h = 2 * R->dpi;
                geo_rect(R, rx - h, ry - h, 2 * h, 2 * h, (SDL_Color){ 110, 110, 125, 255 });
            }
        }
        geo_flush(R);
    }

    /* masses */
    float rad = (float)(sc * 0.17);
    if (rad < 1.5f * R->dpi) rad = 1.5f * R->dpi;
    if (rad > 7 * R->dpi) rad = 7 * R->dpi;
    double vref = 2.0 * sim_omega0(s);
    int defect = sim_defect_index(s);
    if (defect >= 0)
        geo_disc(R, PX(defect), PY(defect), rad * 2.4f,
                 s->p.disorder == DIS_LIGHT ? (SDL_Color){ 120, 255, 160, 200 } : (SDL_Color){ 255, 120, 60, 200 });
    if (s->held >= 0) geo_disc(R, PX(s->held), PY(s->held), rad * 2.2f, (SDL_Color){ 255, 255, 100, 230 });
    for (int i = 0; i < s->n; i++) {
        if (s->fixed[i]) continue;
        double spd = sqrt(s->vx[i] * s->vx[i] + s->vy[i] * s->vy[i]);
        float r = rad;
        if (s->p.disorder == DIS_RANDOM) r = rad * (float)sqrt(s->m[i] / s->p.mass);
        geo_disc(R, PX(i), PY(i), r, speed_colour(spd, vref));
    }
    geo_flush(R);
#undef PX
#undef PY
    SDL_RenderSetClipRect(ren, NULL);
}

/* ───────────────────────────────────────────────────────── spectral panel ── */

static void blit_image(Renderer *R, const float *img, int w, int h, SDL_Rect dst)
{
    if (!img) return;
    if (!R->tex || R->tex_w != w || R->tex_h != h) {
        if (R->tex) SDL_DestroyTexture(R->tex);
        R->tex = SDL_CreateTexture(R->ren, SDL_PIXELFORMAT_ARGB8888, SDL_TEXTUREACCESS_STREAMING, w, h);
        R->tex_w = w; R->tex_h = h;
        if (!R->tex) return;
    }
    void *pix; int pitch;
    if (SDL_LockTexture(R->tex, NULL, &pix, &pitch)) return;
    for (int y = 0; y < h; y++) {
        Uint32 *row = (Uint32 *)((Uint8 *)pix + (size_t)y * pitch);
        for (int x = 0; x < w; x++) row[x] = heat(img[(size_t)y * w + x]);
    }
    SDL_UnlockTexture(R->tex);
    SDL_RenderCopy(R->ren, R->tex, NULL, &dst);
}

static void text(Renderer *R, int x, int y, SDL_Color c, const char *fmt, ...)
{
    char buf[160];
    va_list ap; va_start(ap, fmt); vsnprintf(buf, sizeof buf, fmt, ap); va_end(ap);
    set_col(R->ren, c);
    font_draw(R->ren, x, y, buf, R->fs);
}

static const char *POL_NAME[POL_COUNT] = { "L", "T", "L+T" };

static void plot_hist(Renderer *R, const Hist *h, const int *series, const SDL_Color *cols,
                      int nser, SDL_Rect pr, int include_zero)
{
    if (h->len < 2) return;
    double ymin = include_zero ? 0 : 1e300, ymax = include_zero ? 0 : -1e300;
    for (int i = 0; i < h->len; i++)
        for (int k = 0; k < nser; k++) {
            double v = h->v[(size_t)i * h->nser + series[k]];
            if (v < ymin) ymin = v;
            if (v > ymax) ymax = v;
        }
    if (ymax - ymin < 1e-300) ymax = ymin + 1;
    double t0 = h->t[0], t1 = h->t[h->len - 1];
    if (t1 <= t0) t1 = t0 + 1;
    SDL_FPoint *pts = malloc((size_t)h->len * sizeof *pts);
    if (!pts) return;
    for (int k = 0; k < nser; k++) {
        for (int i = 0; i < h->len; i++) {
            double v = h->v[(size_t)i * h->nser + series[k]];
            pts[i].x = (float)(pr.x + (h->t[i] - t0) / (t1 - t0) * (pr.w - 1));
            pts[i].y = (float)(pr.y + pr.h - 1 - (v - ymin) / (ymax - ymin) * (pr.h - 1));
        }
        set_col(R->ren, cols[k]);
        SDL_RenderDrawLinesF(R->ren, pts, h->len);
    }
    free(pts);
}

static void render_spectrum(Renderer *R, const Sim *s, Spectrum *sp, const UiState *ui)
{
    SDL_Renderer *ren = R->ren;
    SDL_Rect pr = R->spec;
    int lh = font_line_height(R->fs), tx = pr.x, ty = pr.y - lh - 2 * R->fs;
    const SDL_Color title = { 200, 200, 255, 255 }, dim = { 130, 130, 150, 255 };

    SDL_SetRenderDrawColor(ren, 8, 8, 14, 255);
    SDL_RenderFillRect(ren, &pr);

    switch (ui->view) {
    case VIEW_OCC: {
        int w, h; const float *img = spectrum_occ_image(sp, ui->pol, &w, &h);
        blit_image(R, img, w, h, pr);
        text(R, tx, ty, title, "[1/4] |v_%s(q)|^2  occupation", POL_NAME[ui->pol]);
        if (sp->periodic) {
            text(R, pr.x, pr.y + pr.h + 2 * R->fs, dim, "qx -pi..pi, qy -pi..pi (centred)");
        } else {
            text(R, pr.x, pr.y + pr.h + 2 * R->fs, dim, "qx 0..pi ->   qy 0..pi ^   (%s)",
                 sp->bc == BC_CLAMPED ? "DST-I" : "DCT-II");
        }
        break;
    }
    case VIEW_DISP: {
        int w, h; const float *img = spectrum_disp_image(sp, ui->pol, &w, &h);
        blit_image(R, img, w, h, pr);
        double bin = sp->omega_nyq / (sp->nt / 2);
        for (int b = 0; b < 2; b++) {
            if ((b == 0 && ui->pol == POL_T) || (b == 1 && ui->pol == POL_L)) continue;
            SDL_FPoint pts[1100];
            int n = 0;
            for (int ix = 0; ix < sp->nlq && n < 1100; ix++) {
                double wv = b == 0 ? sp->wL[ix] : sp->wT[ix];
                pts[n].x = (float)(pr.x + (ix + 0.5) * pr.w / sp->nlq);
                pts[n].y = (float)(pr.y + pr.h * (h - 0.5 - wv / bin) / h);
                n++;
            }
            SDL_Color c = b == 0 ? COL_L : COL_T; c.a = 230;
            set_col(ren, c);
            for (int i = 0; i + 1 < n; i += 2)        /* dashed */
                SDL_RenderDrawLineF(ren, pts[i].x, pts[i].y, pts[i + 1].x, pts[i + 1].y);
        }
        text(R, tx, ty, title, "[2/4] S_%s(qx,w)  qy=%.2f", POL_NAME[ui->pol], sp->qy[0]);
        int fill = (int)(100.0 * sp->nsamp / sp->nt);
        text(R, pr.x, pr.y + pr.h + 2 * R->fs, dim, "qx 0..pi ->  w 0..%.0f rad/s ^", sp->omega_nyq);
        if (fill < 100) text(R, pr.x + 4 * R->fs, pr.y + 4 * R->fs, (SDL_Color){ 255, 255, 255, 255 },
                             "filling %d%% (%.1f s window)", fill, sp->nt * sp->dt_s);
        text(R, pr.x + pr.w - font_width("-- L  -- T", R->fs) - 4 * R->fs, pr.y + 4 * R->fs, COL_L, "-- L");
        text(R, pr.x + pr.w - font_width("-- T", R->fs) - 4 * R->fs, pr.y + 4 * R->fs, COL_T, "-- T");
        break;
    }
    case VIEW_MODES: {
        int series[2 * SPEC_KMODES]; SDL_Color cols[2 * SPEC_KMODES]; int n = 0;
        for (int b = 0; b < 2; b++) {
            if ((b == 0 && ui->pol == POL_T) || (b == 1 && ui->pol == POL_L)) continue;
            for (int k = 0; k < sp->K; k++) {
                series[n] = b * SPEC_KMODES + k;
                cols[n] = MODE_COL[k];
                if (b == 1 && ui->pol == POL_BOTH) { cols[n].r /= 2; cols[n].g /= 2; cols[n].b /= 2; }
                n++;
            }
        }
        SDL_Rect in = { pr.x + 2, pr.y + 2, pr.w - 4, pr.h - 4 };
        plot_hist(R, &sp->modes, series, cols, n, in, 1);
        text(R, tx, ty, title, "[3/4] mode energy, %s, qy=%.2f", POL_NAME[ui->pol], line_q(sp->bc, s->rows, line_w0(sp->bc)));
        int x = pr.x, y = pr.y + pr.h + 2 * R->fs;
        text(R, x, y, dim, "w="); x += 2 * 6 * R->fs;
        for (int k = 0; k < sp->K; k++) {
            text(R, x, y, MODE_COL[k], "%d", sp->wx0 + k);
            x += 2 * 6 * R->fs;
        }
        double span = sp->modes.len > 1 ? sp->modes.t[sp->modes.len - 1] - sp->modes.t[0] : 0;
        text(R, x + 6 * R->fs, y, dim, "t: %.0f s", span);
        break;
    }
    case VIEW_ENERGY: {
        int series[3] = { 0, 1, 2 };
        SDL_Color cols[3] = { { 90, 160, 255, 255 }, { 255, 110, 110, 255 }, { 235, 235, 235, 255 } };
        SDL_Rect in = { pr.x + 2, pr.y + 2, pr.w - 4, pr.h - 4 };
        plot_hist(R, &sp->energy, series, cols, 3, in, 1);
        text(R, tx, ty, title, "[4/4] energy vs time");
        int y = pr.y + pr.h + 2 * R->fs, x = pr.x;
        text(R, x, y, cols[0], "KE"); x += 4 * 6 * R->fs;
        text(R, x, y, cols[1], "PE"); x += 4 * 6 * R->fs;
        text(R, x, y, cols[2], "E");  x += 3 * 6 * R->fs;
        double span = sp->energy.len > 1 ? sp->energy.t[sp->energy.len - 1] - sp->energy.t[0] : 0;
        text(R, x, y, dim, "t: %.0f s", span);
        break;
    }
    default: break;
    }
    SDL_SetRenderDrawColor(ren, 60, 60, 80, 255);
    SDL_RenderDrawRect(ren, &pr);
}

/* ───────────────────────────────────────────────────────── text panel ── */

static void render_text(Renderer *R, const Sim *s, const UiState *ui)
{
    const int fs = R->fs, lh = font_line_height(fs), x = R->text.x, kx = x + 8 * 6 * fs;
    int y = R->text.y;
    const SDL_Color head = { 200, 200, 255, 255 }, key = { 100, 210, 120, 255 },
                    desc = { 150, 150, 170, 255 }, on = { 80, 200, 100, 255 },
                    off = { 180, 80, 100, 255 }, warn = { 255, 170, 60, 255 };
    const int bottom = R->text.y + R->text.h;

    if (ui->show_help) {
        static const char *help[][2] = {
            { "L-drag", "pull, release to throw" }, { "R-click", "random kick" },
            { "Space",  "pluck centre" },           { "G",      "big random kick" },
            { "R",      "reset" },                  { "B",      "boundary cycle" },
            { "D",      "diagonal springs" },       { "T",      "tension (prestrain)" },
            { "M",      "mass disorder cycle" },    { "A",      "absorbing layer" },
            { "+/-  0", "damping / on-off" },       { "F",      "FPUT preset" },
            { "Sh-F",   "excite lowest mode" },     { "Tab",    "spectrum view" },
            { "V",      "branch L / T / L+T" },     { "P",      "pause" },
            { ", .",    "slower / faster" },        { "H",      "hide this help" },
            { "Q/Esc",  "quit" },
        };
        text(R, x, y, head, "CONTROLS"); y += lh + 2 * fs;
        for (size_t i = 0; i < sizeof help / sizeof *help && y + lh < bottom; i++) {
            text(R, x, y, key, "%s", help[i][0]);
            text(R, kx, y, desc, "%s", help[i][1]);
            y += lh;
        }
        y += lh / 2;
    }
    if (y + lh >= bottom) return;
    text(R, x, y, head, "STATE"); y += lh + 2 * fs;

#define LINE(col, ...) do { if (y + lh < bottom) { text(R, x, y, col, __VA_ARGS__); y += lh; } } while (0)
    LINE(desc, "grid %dx%d  k=%.0f  m=%.2f", s->cols, s->rows, s->p.k, s->p.mass);
    LINE(desc, "boundary  %s", bc_name(s->p.bc));
    LINE(s->p.diagonal ? on : off, "diagonal  %s", s->p.diagonal ? "on" : "off");
    if (s->p.prestrain > 0 && s->prestrain_eff == 0)
        LINE(warn, "tension   %.0f%% (needs clamped/periodic)", 100 * s->p.prestrain);
    else
        LINE(desc, "tension   %.0f%%", 100 * s->prestrain_eff);
    LINE(desc, "damping   %.3g kg/s", s->p.damping);
    LINE(s->p.absorb ? on : desc, "absorber  %s", s->p.absorb ? "on" : "off");
    LINE(desc, "masses    %s", disorder_name(s->p.disorder));
    LINE(desc, "speed     %gx%s", ui->speed, ui->paused ? "  PAUSED" : "");
    LINE(desc, "t = %.1f s   %.0f fps", s->t, ui->fps);
    LINE(desc, "E = %.4g", ui->E);
    if (ui->drift_ok && fabs(ui->E0) > 1e-12)
        LINE(desc, "dE/E = %+.2e", (ui->E - ui->E0) / ui->E0);
    else
        LINE(desc, "dE/E   n/a (damped/driven)");
    double wdt = s->omega_max * s->p.dt;
    if (wdt > 0.5) LINE(warn, "w_max*dt = %.2f: reduce --dt", wdt);
    if (ui->lagging) LINE(warn, "CPU-bound: running < %gx", ui->speed);
#undef LINE
}

/* ───────────────────────────────────────────────────────── frame ── */

void render_frame(Renderer *R, const Sim *s, Spectrum *sp, const UiState *ui)
{
    SDL_Renderer *ren = R->ren;
    SDL_SetRenderDrawColor(ren, 15, 15, 20, 255);
    SDL_RenderClear(ren);

    render_lattice(R, s);

    SDL_SetRenderDrawColor(ren, 12, 12, 18, 255);
    SDL_Rect panel = { R->sim.w, 0, R->out_w - R->sim.w, R->out_h };
    SDL_RenderFillRect(ren, &panel);
    SDL_SetRenderDrawColor(ren, 40, 40, 60, 255);
    SDL_RenderDrawLine(ren, R->sim.w, 0, R->sim.w, R->out_h);

    render_spectrum(R, s, sp, ui);
    render_text(R, s, ui);

    if (ui->now < ui->toast_until && ui->toast[0]) {
        int w = font_width(ui->toast, R->fs), lh = font_line_height(R->fs);
        SDL_Rect bg = { R->sim.x + 8 * R->fs, R->sim.h - lh - 12 * R->fs, w + 8 * R->fs, lh + 4 * R->fs };
        SDL_SetRenderDrawColor(ren, 0, 0, 0, 190);
        SDL_RenderFillRect(ren, &bg);
        text(R, bg.x + 4 * R->fs, bg.y + 3 * R->fs, (SDL_Color){ 255, 230, 140, 255 }, "%s", ui->toast);
    }
}
