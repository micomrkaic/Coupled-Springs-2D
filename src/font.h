/*
 * font.h — 5x7 bitmap font (ASCII 32..126), batched into SDL_RenderFillRects
 */
#ifndef FONT_H
#define FONT_H

#include <SDL.h>

/* draws in the renderer's current colour; returns the advance in pixels */
int  font_draw(SDL_Renderer *ren, int x, int y, const char *s, int scale);
int  font_width(const char *s, int scale);
static inline int font_line_height(int scale) { return 10 * scale; }

#endif
