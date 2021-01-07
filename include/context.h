#ifndef CONTEXT_H
#define CONTEXT_H

#include "linmath_d.h"
#include <pthread.h>

typedef unsigned char color_t[3];

typedef struct
{
    size_t size;
    vec3 centroid;
    double (*verts)[3];
    double (*texcoords)[2];
    double (*normals)[3];

}Mesh;

typedef struct
{
    int rows;
    int cols;
    int grid_rows;
    int grid_cols;
    int n_threads;
    int *thread_sizes;
    color_t (*buffer);
    color_t (*buffer_2);  // one buffer for drawing, one for display
    double *z_buffer;
    double *z_buffer_2;
    vec2 **grid_points;
    mat4x4 view_mat;
    Mesh *mesh;

}RenderCtx;

pthread_t clear_buffers_start(RenderCtx *ctx);
void wait_for_clear(pthread_t thread);
void swap_buffers(RenderCtx *ctx);
void move_camera(RenderCtx *ctx, vec3 offset);

RenderCtx init_renderer();
void destroy_renderer(RenderCtx *ctx);

#endif
