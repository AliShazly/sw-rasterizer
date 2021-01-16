#ifndef DRAW_H
#define DRAW_H

#include "utils.h"
#include "linmath_d.h"
#include "context.h"

void draw_triangle(vec4 triangle[TRI_NPTS], color_t fill, bool two_sided, RenderCtx *ctx);
void draw_line(vec2 p1, vec2 p2, color_t color, int thick, RenderCtx *ctx);
void draw_point(vec2 p, color_t color, int radius, RenderCtx *ctx);

#endif
