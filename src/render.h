/*
 * render.h — SDL2 drawing: lattice, spectral panel, help/state text
 *
 * All drawing is in renderer output pixels, so HiDPI displays get full
 * resolution; window (point) coordinates from mouse events are converted
 * with render_to_lattice().  Springs and masses are batched into one
 * SDL_RenderGeometry call each (needs SDL >= 2.0.18).
 */
#ifndef RENDER_H
#define RENDER_H

#include "sim.h"
#include "spectrum.h"

#include <SDL.h>

typedef struct {
    View        view;
    Pol         pol;
    int         paused, lagging, show_help;
    double      speed;
    double      E, E0;          /* current and post-event energy */
    int         drift_ok;       /* conservative and undisturbed since E0 */
    double      fps;
    char        toast[96];
    double      toast_until, now;
} UiState;

typedef struct {
    SDL_Window   *win;
    SDL_Renderer *ren;
    int           out_w, out_h, win_w, win_h;
    float         dpi;
    int           fs;              /* font scale */
    SDL_Rect      sim, spec, text;
    double        scale, ox, oy;   /* pixel = o + lattice * scale */

    SDL_Texture  *tex;
    int           tex_w, tex_h;

    SDL_Vertex   *vtx;
    int          *idx;
    int           vcap, icap, nv, ni;
} Renderer;

int  render_init(Renderer *R, SDL_Window *win, SDL_Renderer *ren);
void render_free(Renderer *R);
void render_layout(Renderer *R, const Sim *s);
/* draws the frame; the caller presents it */
void render_frame(Renderer *R, const Sim *s, Spectrum *sp, const UiState *ui);

int  render_in_sim(const Renderer *R, int wx, int wy);
void render_to_lattice(const Renderer *R, int wx, int wy, double *lx, double *ly);

#endif
