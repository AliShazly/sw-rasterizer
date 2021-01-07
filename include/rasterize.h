#ifndef RASTERIZE_H
#define RASTERIZE_H

#include "context.h"

void draw_object_threads(RenderCtx *ctx);
void draw_object_wireframe(RenderCtx *ctx);
void draw_grid(RenderCtx *ctx);

// used for precompution
vec2 **compute_grid(double size, int n_subdivs, int *out_rows, int *out_cols);
void lookat(vec3 camera_pos, vec3 target, mat4x4 dst, vec3 up);

#endif
