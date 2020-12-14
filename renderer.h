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
    vec3 camera_pos;
    vec3 camera_z;
    mat4x4 cam_space_transform;
    uint8_t (*buffer)[3];
    double *z_buffer;
    uint8_t *new_z;
    ObjMesh *mesh;

}RenderCtx;


double inc(double i, double lim);
int clamp(int val, int min, int max);
double normalize(double val, double upper, double lower);
void normalize_vec3(vec3 dst, vec3 val, double upper, double lower);
double gamma_correct(double val, double g);
double distance(vec3 a, vec3 b);

void mesh_bounds(vec3 *verts, size_t n_verts, vec3 out_min, vec3 out_max);
void mesh_centroid(vec3 dst, vec3 *verts, size_t n_verts);
void normalize_coords(vec3 *out_arr, vec3 *verts, size_t n_verts);

int orient2d(vec2 a, vec2 b, vec2 c);
void triangle_bbox(vec3 triangle[3], double *max_x, double *max_y, double *min_x, double *min_y);
void draw_triangle(vec3 triangle[3], uint8_t fill[3], RenderCtx *ctx);

void surface_normal(vec3 dst, vec3 *triangle);
void camera_transform(vec3 camera_pos, vec3 target, vec3 up, mat4x4 dst, vec3 out_z);
void ortho_projection(vec2 dst, vec3 pt, vec2 scale, vec2 offset);
void draw_object(double (*verts)[3], size_t n_verts, RenderCtx *ctx);

RenderCtx init_renderer();
void destroy_renderer(RenderCtx *ctx);

#endif

