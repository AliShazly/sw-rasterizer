#ifndef CONTEXT_H
#define CONTEXT_H

#include "linmath_d.h"

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
    vec3 camera_pos;
    vec3 camera_z;
    mat4x4 cam_transform;
    color_t (*buffer);
    double *z_buffer;
    vec2 **grid_points;
    Mesh *mesh;

}RenderCtx;

RenderCtx init_renderer();
void destroy_renderer(RenderCtx *ctx);

#endif
