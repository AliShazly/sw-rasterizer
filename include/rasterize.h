#ifndef RASTERIZE_H
#define RASTERIZE_H

#include "context.h"

void clear_buffers(RenderCtx *ctx);
void draw_object_threads(RenderCtx *ctx);
void draw_object_wireframe(RenderCtx *ctx);
void draw_grid(RenderCtx *ctx);

// used for precompution
void camera_transform(vec3 camera_pos, vec3 target, vec3 up, mat4x4 dst, vec3 out_z);
vec2 **compute_grid(double size, int n_subdivs, int *out_rows, int *out_cols);

#endif
