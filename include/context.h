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
    color_t (*buffer);
    double *z_buffer;
    vec2 **grid_points;
    mat4x4 view_mat;
    Mesh *mesh;

}RenderCtx;

void clear_buffers(RenderCtx *ctx);
void move_camera(RenderCtx *ctx, vec3 offset);

RenderCtx init_renderer();
void destroy_renderer(RenderCtx *ctx);

#endif
