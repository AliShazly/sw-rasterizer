#ifndef RENDERER_H
#define RENDERER_H

#include <stddef.h>
#include <stdint.h>

#include "lib/linmath_d.h"

typedef struct
{
    size_t size;
    vec3 centroid;
    double (*verts)[3];
    double (*texcoords)[2];
    double (*normals)[3];

}ObjMesh;

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
    uint8_t (*buffer)[3];
    double *z_buffer;
    vec2 **grid_points;
    ObjMesh *mesh;

}RenderCtx;

typedef struct
{
    size_t start_idx;
    size_t end_idx;
    double time; // temporary
    RenderCtx *ctx;

}ThreadArgs;


double inc(double i, double lim);
int clamp(int val, int min, int max);
double normalize(double val, double upper, double lower);
double gamma_correct(double val, double g);
double distance(vec3 a, vec3 b);
void line_midpoint(vec2 dst, vec2 a, vec2 b);
void subdivide_line(vec2 *dst, size_t out_len, vec2 line[2], int n_subdivs);

void mesh_bounds(vec3 *verts, size_t n_verts, vec3 out_min, vec3 out_max);
void mesh_centroid(vec3 dst, vec3 *verts, size_t n_verts);
void surface_normal(vec3 dst, vec3 *triangle);
void divide_among_threads(int data_size, int n_threads, int out_sizes[n_threads]);
double light_triangle(vec3 camera_z, vec3 world_tri[3]);
void normalize_coords(vec3 *out_arr, vec3 *verts, size_t n_verts);

int orient2d(vec2 a, vec2 b, vec2 c);
void triangle_bbox(vec3 triangle[3], double *max_x, double *max_y, double *min_x, double *min_y);
void draw_triangle(vec3 triangle[3], uint8_t fill[3], RenderCtx *ctx);

void camera_transform(vec3 camera_pos, vec3 target, vec3 up, mat4x4 dst, vec3 out_z);
void ortho_projection(vec2 dst, vec3 pt, vec2 scale, vec2 offset);
void world_to_screen(vec3 dst, mat4x4 cam_space_transform, vec4 world_pt, int rows, int cols);
void clear_buffers(RenderCtx *ctx);
void *object_thread(void *thread_args);
void draw_object_threads(RenderCtx *ctx);
void draw_object_wireframe(RenderCtx *ctx);
vec2 **compute_grid(double size, int n_subdivs, int *out_rows, int *out_cols);
void draw_grid(RenderCtx *ctx);

RenderCtx init_renderer();
void destroy_renderer(RenderCtx *ctx);

#endif

