/*
 * spring_grid — real-time 2-D spring-mass lattice with spectral diagnostics
 *
 * Window:  [ lattice | spectrum view + controls/state ]
 * Physics: sim.c (lattice units, split-damped velocity Verlet)
 * Spectra: spectrum.c (FFTW), analytic branches: modes.c
 * Controls: see the help panel or README.md;  options: spring_grid --help
 */
#define _POSIX_C_SOURCE 200809L   /* strdup, strtok_r */

#include "params.h"
#include "render.h"
#include "sim.h"
#include "spectrum.h"
#include "util.h"

#include <SDL.h>

#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
    Sim       sim;
    Spectrum  sp;
    Renderer  R;
    UiState   ui;
    MotionLog mlog;
    double    saved_damping;
    int       running;
} App;

static double now_s(void)
{
    return (double)SDL_GetPerformanceCounter() / (double)SDL_GetPerformanceFrequency();
}

static void toast(App *a, const char *fmt, ...)
{
    va_list ap; va_start(ap, fmt);
    vsnprintf(a->ui.toast, sizeof a->ui.toast, fmt, ap);
    va_end(ap);
    a->ui.toast_until = now_s() + 2.5;
}

/* energy baseline for the drift readout: reset after anything that adds or
 * removes energy from outside */
static void mark_event(App *a)
{
    double ke, pe; sim_energy(&a->sim, &ke, &pe);
    a->ui.E0 = ke + pe;
}

static void structural_change(App *a)
{
    sim_rebuild(&a->sim);
    if (spectrum_rebuild(&a->sp, &a->sim)) { fprintf(stderr, "spectrum: out of memory\n"); exit(1); }
    mark_event(a);
}

static int can_periodic(const Sim *s) { return s->cols % 2 == 0 && s->rows % 2 == 0; }

/* ───────────────────────────────────────────────────────── keys ── */

static void handle_key(App *a, SDL_Keycode key, Uint16 mod)
{
    Sim *s = &a->sim;
    Params *p = &s->p;
    UiState *ui = &a->ui;
    const int shift = (mod & KMOD_SHIFT) != 0;
    const double w0 = sim_omega0(s);

    switch (key) {
    case SDLK_ESCAPE: case SDLK_q: a->running = 0; break;

    case SDLK_r:
        sim_reset(s); spectrum_clear_history(&a->sp, s); mark_event(a);
        break;
    case SDLK_SPACE:
        sim_displace(s, sim_idx(s, s->cols / 2, s->rows / 2), 0.0, -2.0);
        mark_event(a);
        break;
    case SDLK_g:
        sim_kick_random(s, 0.5 * w0, 0.5 * w0);
        mark_event(a);
        break;

    case SDLK_b: {
        Boundary bc = p->bc;
        do bc = (Boundary)((bc + 1) % BC_COUNT); while (bc == BC_PERIODIC && !can_periodic(s));
        p->bc = bc;
        structural_change(a);
        toast(a, "boundary: %s%s", bc_name(bc),
              bc == BC_PINNED ? " (edges tied to rest sites)" : bc == BC_CLAMPED ? " (edges held fixed)" : "");
        break;
    }
    case SDLK_d:
        p->diagonal = !p->diagonal;
        structural_change(a);
        toast(a, "diagonal springs %s", p->diagonal ? "on" : "off");
        break;
    case SDLK_t: {
        static const double steps[] = { 0.0, 0.02, 0.05, 0.10 };
        int k = 0;
        while (k < 4 && steps[k] <= p->prestrain + 1e-9) k++;
        p->prestrain = steps[k % 4];
        structural_change(a);
        if (p->prestrain > 0 && s->prestrain_eff == 0)
            toast(a, "prestrain %.0f%% set; acts only under clamped/periodic (B)", 100 * p->prestrain);
        else
            toast(a, "prestrain %.0f%%", 100 * p->prestrain);
        break;
    }
    case SDLK_m:
        p->disorder = (Disorder)((p->disorder + 1) % DIS_COUNT);
        sim_set_masses(s);
        structural_change(a);
        toast(a, "masses: %s", p->disorder == DIS_LIGHT ? "light defect at centre"
                             : p->disorder == DIS_HEAVY ? "heavy defect at centre"
                             : p->disorder == DIS_RANDOM ? "random" : "uniform");
        break;
    case SDLK_a:
        p->absorb = !p->absorb;
        sim_update_damping(s);
        mark_event(a);
        toast(a, "absorbing layer %s (%d nodes)", p->absorb ? "on" : "off", p->absorb_width);
        break;

    case SDLK_PLUS: case SDLK_EQUALS: case SDLK_KP_PLUS:
        p->damping = p->damping > 0 ? p->damping * 1.25 : 0.01;
        sim_update_damping(s); mark_event(a);
        break;
    case SDLK_MINUS: case SDLK_KP_MINUS:
        p->damping /= 1.25;
        if (p->damping < 1e-3) p->damping = 0.0;
        sim_update_damping(s); mark_event(a);
        break;
    case SDLK_0: case SDLK_KP_0:
        if (p->damping > 0) { a->saved_damping = p->damping; p->damping = 0.0; }
        else p->damping = a->saved_damping > 0 ? a->saved_damping : 0.4;
        sim_update_damping(s); mark_event(a);
        toast(a, "damping %.3g", p->damping);
        break;

    case SDLK_f:
        if (shift) {                        /* lowest mode, current settings */
            int b = ui->pol == POL_T ? 1 : 0;
            sim_excite_mode(s, b, p->fput_strain);
            spectrum_clear_history(&a->sp, s);
            ui->view = VIEW_MODES;
            mark_event(a);
            toast(a, "lowest %s mode, strain %.2f", b ? "T" : "L", p->fput_strain);
        } else {                            /* FPUT preset */
            if (!can_periodic(s)) { toast(a, "FPUT preset needs even cols and rows"); break; }
            p->bc = BC_PERIODIC; p->diagonal = 1; p->prestrain = 0.0;
            if (p->damping > 0) a->saved_damping = p->damping;
            p->damping = 0.0; p->absorb = 0;
            if (p->disorder != DIS_NONE) { p->disorder = DIS_NONE; sim_set_masses(s); }
            sim_rebuild(s);
            sim_excite_mode(s, 0, p->fput_strain);
            s->t = 0; s->steps = 0;
            if (spectrum_rebuild(&a->sp, s)) exit(1);
            mark_event(a);
            ui->view = VIEW_MODES; ui->pol = POL_L; ui->speed = 4.0;
            toast(a, "FPUT: periodic, diagonals, L mode w=1, strain %.2f, 4x", p->fput_strain);
        }
        break;

    case SDLK_TAB:
        ui->view = (View)((ui->view + (shift ? VIEW_COUNT - 1 : 1)) % VIEW_COUNT);
        break;
    case SDLK_v:
        ui->pol = (Pol)((ui->pol + 1) % POL_COUNT);
        break;
    case SDLK_p:
        ui->paused = !ui->paused;
        break;
    case SDLK_COMMA:
        if (ui->speed > 1.0 / 16) ui->speed /= 2;
        toast(a, "speed %gx", ui->speed);
        break;
    case SDLK_PERIOD:
        if (ui->speed < 64) ui->speed *= 2;
        toast(a, "speed %gx", ui->speed);
        break;
    case SDLK_h:
        ui->show_help = !ui->show_help;
        break;
    default: break;
    }
}

/* ───────────────────────────────────────────────────────── mouse ── */

static void drag_target(App *a, int wx, int wy, double t)
{
    Sim *s = &a->sim;
    double lx, ly;
    render_to_lattice(&a->R, wx, wy, &lx, &ly);
    mlog_push(&a->mlog, t, lx, ly);
    s->held_ux = lx - s->held % s->cols;
    s->held_uy = ly - s->held / s->cols;
}

static void handle_mouse(App *a, const SDL_Event *ev)
{
    Sim *s = &a->sim;
    double t = now_s();
    switch (ev->type) {
    case SDL_MOUSEBUTTONDOWN:
        if (!render_in_sim(&a->R, ev->button.x, ev->button.y)) break;
        if (ev->button.button == SDL_BUTTON_LEFT) {
            double lx, ly;
            render_to_lattice(&a->R, ev->button.x, ev->button.y, &lx, &ly);
            int i = sim_nearest(s, lx, ly, 0.9);
            if (i >= 0) {
                s->held = i;
                s->held_vx = s->held_vy = 0.0;
                mlog_clear(&a->mlog);
                drag_target(a, ev->button.x, ev->button.y, t);
            }
        } else if (ev->button.button == SDL_BUTTON_RIGHT) {
            double w0 = sim_omega0(s);
            sim_kick_random(s, 0.25 * w0, 0.75 * w0);
            mark_event(a);
        }
        break;
    case SDL_MOUSEMOTION:
        if (s->held >= 0) drag_target(a, ev->motion.x, ev->motion.y, t);
        break;
    case SDL_MOUSEBUTTONUP:
        if (ev->button.button == SDL_BUTTON_LEFT && s->held >= 0) {
            drag_target(a, ev->button.x, ev->button.y, t);
            double vx, vy;
            mlog_velocity(&a->mlog, t, 0.08, &vx, &vy);
            int i = s->held;
            s->held = -1;
            s->vx[i] = vx; s->vy[i] = vy;
            mark_event(a);
        }
        break;
    }
}

/* ───────────────────────────────────────────────────────── scripted input ──
 * SG_SCRIPT="frame:action[:arg];..."   actions: key:<SDL key name>[+shift],
 * mdown:x,y  mmove:x,y  mup:x,y (left button, window points), shot:<file.bmp>,
 * quit.  Used for headless screenshots and smoke tests. */

typedef struct { int frame; char act[8]; char arg[200]; } ScriptItem;

static int load_script(ScriptItem *items, int cap)
{
    const char *env = getenv("SG_SCRIPT");
    if (!env) return 0;
    char *buf = strdup(env), *save = NULL;
    int n = 0;
    for (char *tok = strtok_r(buf, ";", &save); tok && n < cap; tok = strtok_r(NULL, ";", &save)) {
        ScriptItem it = { 0 };
        char *c1 = strchr(tok, ':');
        if (!c1) continue;
        *c1 = 0; it.frame = atoi(tok);
        char *c2 = strchr(c1 + 1, ':');
        if (c2) { *c2 = 0; snprintf(it.arg, sizeof it.arg, "%s", c2 + 1); }
        snprintf(it.act, sizeof it.act, "%s", c1 + 1);
        items[n++] = it;
    }
    free(buf);
    return n;
}

static void save_shot(App *a, const char *path)
{
    SDL_Surface *surf = SDL_CreateRGBSurfaceWithFormat(0, a->R.out_w, a->R.out_h, 32, SDL_PIXELFORMAT_ARGB8888);
    if (!surf) return;
    if (SDL_RenderReadPixels(a->R.ren, NULL, SDL_PIXELFORMAT_ARGB8888, surf->pixels, surf->pitch) == 0)
        SDL_SaveBMP(surf, path);
    fprintf(stderr, "shot %s at sim t = %.2f s\n", path, a->sim.t);
    SDL_FreeSurface(surf);
}

/* ───────────────────────────────────────────────────────── main ── */

int main(int argc, char **argv)
{
    Params p;
    params_default(&p);
    int pr = params_parse(&p, argc, argv);
    if (pr) return pr > 0 ? 0 : 2;

    if (SDL_Init(SDL_INIT_VIDEO)) { fprintf(stderr, "SDL_Init: %s\n", SDL_GetError()); return 1; }
    SDL_SetHint(SDL_HINT_RENDER_SCALE_QUALITY, "0");

    int ww = p.win_w ? p.win_w : 1200, wh = p.win_h ? p.win_h : 800;
    SDL_Rect ub;
    if (!p.win_w && !p.win_h && SDL_GetDisplayUsableBounds(0, &ub) == 0) {
        if (ww > ub.w * 92 / 100) ww = ub.w * 92 / 100;
        if (wh > ub.h * 92 / 100) wh = ub.h * 92 / 100;
    }
    SDL_Window *win = SDL_CreateWindow("Spring-Mass Grid", SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED,
                                       ww, wh, SDL_WINDOW_RESIZABLE | SDL_WINDOW_ALLOW_HIGHDPI);
    if (!win) { fprintf(stderr, "SDL_CreateWindow: %s\n", SDL_GetError()); return 1; }
    SDL_SetWindowMinimumSize(win, 640, 480);
    SDL_Renderer *ren = SDL_CreateRenderer(win, -1, SDL_RENDERER_ACCELERATED | SDL_RENDERER_PRESENTVSYNC);
    if (!ren) ren = SDL_CreateRenderer(win, -1, SDL_RENDERER_SOFTWARE);
    if (!ren) { fprintf(stderr, "SDL_CreateRenderer: %s\n", SDL_GetError()); return 1; }
    SDL_SetRenderDrawBlendMode(ren, SDL_BLENDMODE_BLEND);

    static App app;
    App *a = &app;
    if (sim_init(&a->sim, &p) || spectrum_init(&a->sp, &a->sim) || render_init(&a->R, win, ren)) {
        fprintf(stderr, "out of memory\n"); return 1;
    }
    if (a->sim.omega_max * p.dt > 0.5)
        fprintf(stderr, "warning: omega_max*dt = %.2f; the integrator is inaccurate above ~0.5 "
                        "and unstable at 2. Reduce --dt.\n", a->sim.omega_max * p.dt);
    render_layout(&a->R, &a->sim);
    a->ui = (UiState){ .view = VIEW_OCC, .pol = POL_BOTH, .speed = p.speed, .show_help = 1 };
    a->saved_damping = p.damping > 0 ? p.damping : 0.4;
    a->running = 1;

    sim_displace(&a->sim, sim_idx(&a->sim, a->sim.cols / 2, a->sim.rows / 2), 0.0, -2.5);
    mark_event(a);

    ScriptItem script[64];
    int nscript = load_script(script, 64), si = 0;
    long frame = 0;

    double t_prev = now_s(), acc = 0, fps_t0 = t_prev;
    double lag_wall0 = t_prev, lag_sim0 = 0, lag_want = 0;   /* 1 s throughput window */
    int fps_n = 0;
    SDL_Event ev;

    while (a->running) {
        while (SDL_PollEvent(&ev)) {
            switch (ev.type) {
            case SDL_QUIT: a->running = 0; break;
            case SDL_KEYDOWN: handle_key(a, ev.key.keysym.sym, ev.key.keysym.mod); break;
            case SDL_MOUSEBUTTONDOWN: case SDL_MOUSEBUTTONUP: case SDL_MOUSEMOTION:
                handle_mouse(a, &ev); break;
            case SDL_WINDOWEVENT:
                if (ev.window.event == SDL_WINDOWEVENT_SIZE_CHANGED) render_layout(&a->R, &a->sim);
                break;
            }
        }
        while (si < nscript && script[si].frame <= frame) {
            ScriptItem *it = &script[si++];
            if (!strcmp(it->act, "key")) {
                char name[200]; snprintf(name, sizeof name, "%s", it->arg);
                Uint16 mod = 0; char *plus = strstr(name, "+shift");
                if (plus) { *plus = 0; mod = KMOD_LSHIFT; }
                handle_key(a, SDL_GetKeyFromName(name), mod);
            } else if (it->act[0] == 'm') {
                SDL_Event e = { 0 };
                int x = 0, y = 0;
                sscanf(it->arg, "%d,%d", &x, &y);
                if (!strcmp(it->act, "mmove")) {
                    e.type = SDL_MOUSEMOTION; e.motion.x = x; e.motion.y = y;
                } else {
                    e.type = strcmp(it->act, "mdown") ? SDL_MOUSEBUTTONUP : SDL_MOUSEBUTTONDOWN;
                    e.button.button = SDL_BUTTON_LEFT; e.button.x = x; e.button.y = y;
                }
                handle_mouse(a, &e);
            } else if (!strcmp(it->act, "quit")) {
                a->running = 0;
            }
        }

        /* fixed-timestep accumulator: simulated time advances at `speed`
         * times wall time, independent of the display refresh rate */
        double t = now_s(), fdt = t - t_prev;
        t_prev = t;
        if (fdt > 0.1) fdt = 0.1;
        Sim *s = &a->sim;
        if (!a->ui.paused) acc += fdt * a->ui.speed;
        if (!a->ui.paused) lag_want += fdt * a->ui.speed;
        /* at most 0.1 s of wall time's worth of steps per frame; if the CPU
         * cannot keep up, the run slows down instead of spiralling */
        int max_steps = (int)ceil(a->ui.speed * 0.1 / s->p.dt) + 1, n = 0;

        if (s->held >= 0) {
            double vx, vy;
            mlog_velocity(&a->mlog, t, 0.08, &vx, &vy);
            s->held_vx = vx; s->held_vy = vy;
        }
        while (acc >= s->p.dt && n < max_steps) {
            sim_step(s);
            spectrum_after_step(&a->sp, s);
            acc -= s->p.dt;
            n++;
        }
        if (acc >= s->p.dt) acc = 0;
        if (n) spectrum_frame(&a->sp, s);

        double ke, pe;
        sim_energy(s, &ke, &pe);
        a->ui.E = ke + pe;
        a->ui.drift_ok = s->p.damping == 0 && !s->p.absorb && s->held < 0;
        if (t - lag_wall0 >= 1.0) {
            double got = s->t - lag_sim0;
            a->ui.lagging = lag_want > 0 && got >= 0 && got < 0.9 * lag_want;
            lag_wall0 = t; lag_sim0 = s->t; lag_want = 0;
        }
        a->ui.now = t;
        if (++fps_n, t - fps_t0 > 0.5) { a->ui.fps = fps_n / (t - fps_t0); fps_n = 0; fps_t0 = t; }

        render_frame(&a->R, s, &a->sp, &a->ui);
        for (int k = si - 1; k >= 0 && script[k].frame == frame; k--)
            if (!strcmp(script[k].act, "shot")) save_shot(a, script[k].arg);
        SDL_RenderPresent(ren);
        frame++;
    }

    spectrum_free(&a->sp);
    sim_free(&a->sim);
    render_free(&a->R);
    SDL_DestroyRenderer(ren);
    SDL_DestroyWindow(win);
    SDL_Quit();
    fftw_cleanup();
    return 0;
}
