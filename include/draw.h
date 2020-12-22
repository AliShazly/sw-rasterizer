#ifndef DRAW_H
#define DRAW_H

#include "linmath_d.h"
#include "context.h"

void draw_triangle(vec3 triangle[3], color_t fill, RenderCtx *ctx);
void draw_line(int x0, int y0, int x1, int y1, color_t color, RenderCtx *ctx);

#endif
