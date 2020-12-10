#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <stdint.h>
#include <assert.h>
#include <limits.h>
#include <math.h>
#include <string.h>
#include <time.h>

// https://github.com/datenwolf/linmath.h
//   i made all the vec# types doubles, added string.h
#include "lib/linmath_d.h"

#include "obj_parser.h"
#include "renderer.h"

double normalize(double val, double upper, double lower)
{
    return (val - lower) / (upper - lower);
}

// find bounding box of verts and equalizes the width / height
void mesh_bounds(vec3 *verts, size_t n_verts, vec3 out_min, vec3 out_max)
{
    double min_x, min_y, min_z;
    double max_x, max_y, max_z;
    min_x = min_y = min_z = INT_MAX;
    max_x = max_y = max_z = INT_MIN;

    for (int i = 0; i < n_verts; i++)
    {
        if (verts[i][0] < min_x)
        {
            min_x = verts[i][0];
        }
        else if (verts[i][0] > max_x)
        {
            max_x = verts[i][0];
        }

        if (verts[i][1] < min_y)
        {
            min_y = verts[i][1];
        }
        else if (verts[i][1] > max_y)
        {
            max_y = verts[i][1];
        }

        if (verts[i][2] < min_z)
        {
            min_z = verts[i][2];
        }
        else if (verts[i][2] > max_z)
        {
            max_z = verts[i][2];
        }
    }

    // everything should have a value
    assert(min_x != INT_MAX && min_y != INT_MAX && min_z != INT_MAX);
    assert(max_x != INT_MIN && max_y != INT_MIN && max_z != INT_MIN);

    double width = max_x - min_x;
    double height = max_y - min_y;

    // Extending the smaller two dimensions to form a square x/y plane
    double largest = fmax(width, height);
    max_x += largest - width;
    max_y += largest - height;

    vec3 mins = {min_x, min_y, min_z};
    vec3 maxs = {max_x, max_y, max_z};
    memcpy(out_min, mins, sizeof(vec3));
    memcpy(out_max, maxs, sizeof(vec3));
}

void mesh_centroid(vec3 dst, vec3 *verts, size_t n_verts)
{
    long sum_x = 0, sum_y = 0, sum_z = 0;
    for (int i = 0; i < n_verts; i++)
    {
        sum_x += verts[i][0];
        sum_y += verts[i][1];
        sum_z += verts[i][2];
    }
    dst[0] = sum_x / (double)n_verts;
    dst[1] = sum_y / (double)n_verts;
    dst[2] = sum_z / (double)n_verts;
}

int orient2d(vec2 a, vec2 b, vec2 c)
{
    return (b[0]-a[0])*(c[1]-a[1]) - (b[1]-a[1])*(c[0]-a[0]);
}

// should actually take a vec2 triangle[3] but vec3 is used here because of the packed Z buffer value
void triangle_bbox(vec3 triangle[3], double *max_x, double *max_y, double *min_x, double *min_y)
{
    double *a = triangle[0];
    double *b = triangle[1];
    double *c = triangle[2];

    *max_x = fmax(a[0], fmax(b[0], c[0]));
    *max_y = fmax(a[1], fmax(b[1], c[1]));
    *min_x = fmin(a[0], fmin(b[0], c[0]));
    *min_y = fmin(a[1], fmin(b[1], c[1]));
}

// find the bbox of a triangle, loop over each pixel in bbox,
//    see if it's in the triangle by checking for negative bary coords
// TODO: this only works for square dimensions, and switching the VLA dimensions fixes it??
void draw_triangle(vec3 triangle[3], uint8_t fill[3], RenderCtx *ctx)
{
    double *a = triangle[0];
    double *b = triangle[1];
    double *c = triangle[2];

    double max_x, max_y, min_x, min_y;
    triangle_bbox(triangle, &max_x, &max_y, &min_x, &min_y);

    vec3 sp; // screen space point
    // looping through every point in the bounding box with subpixel accuracy to make sure triangles line up
    // TODO: would only need to loop once per pixel (2x speedup) if this used a 'fill convention'
    for (sp[0] = min_x; sp[0] <= max_x; sp[0] += ctx->subpixel)
    {
        for (sp[1] = min_y; sp[1] <= max_y; sp[1] += ctx->subpixel)
        {
            // if bbox point is out of screen bounds
            if (sp[0] >= ctx->cols - 1 || sp[1] >= ctx->rows - 1 || sp[0] < 0 || sp[1] < 0)
            {
                continue;
            }

            // barycentric coords
            int w0 = orient2d(b, c, sp);
            int w1 = orient2d(c, a, sp);
            int w2 = orient2d(a, b, sp);

            // if point is inside or on all edges
            if (w0 >= 0 && w1 >= 0 && w2 >= 0)
            {
                // init value for Z buffer
                sp[2] = 0;

                // smearing actual z values over bayercentric coords
                sp[2] += a[2] * w0;
                sp[2] += b[2] * w1;
                sp[2] += c[2] * w2;

                // flipping vertically
                int row = (ctx->rows - 1) - ceil(sp[1]);
                int col = ceil(sp[0]);

                // 2D index (row,col) to a 1D index
                int idx = row * ctx->cols + col;

                if (ctx->z_buffer[idx] < sp[2])
                {
                    ctx->z_buffer[idx] = sp[2];
                    ctx->buffer[idx][0] = fill[0];
                    ctx->buffer[idx][1] = fill[1];
                    ctx->buffer[idx][2] = fill[2];
                }
            }
        }
    }
}

// computes surface normal of triangle
void surface_normal(vec3 dst, vec3 *triangle)
{
    vec3 v;
    vec3 w;
    vec3_sub(v, triangle[2], triangle[0]);
    vec3_sub(w, triangle[1], triangle[0]);
    vec3_mul_cross(dst, v, w);
}

// https://www.3dgep.com/understanding-the-view-matrix/#Transformations
void camera_transform(mat4x4 dst, vec3 camera_pos, vec3 target, vec3 up)
{
    vec3 z_axis; // camera forward vector
    vec3_sub(z_axis, camera_pos, target);
    vec3_norm(z_axis, z_axis);

    vec3 x_axis; // camera right vector
    vec3_mul_cross(x_axis, up, z_axis);
    vec3_norm(x_axis, x_axis);

    vec3 y_axis; // camera up vector
    vec3_mul_cross(y_axis, z_axis, x_axis);

    mat4x4 orientation = {
        {x_axis[0], y_axis[0], z_axis[0], 0},
        {x_axis[1], y_axis[1], z_axis[1], 0},
        {x_axis[2], y_axis[2], z_axis[2], 0},
        {0,         0,         0,         1}
    };

    mat4x4 translation;
    mat4x4_identity(translation);
    translation[3][0] = -camera_pos[0];
    translation[3][1] = -camera_pos[1];
    translation[3][2] = -camera_pos[2];

    mat4x4_mul(dst, orientation, translation);
}

// https://en.wikipedia.org/wiki/3D_projection#Orthographic_projection
void ortho_projection(vec2 dst, vec3 pt, vec2 scale, vec2 offset)
{
    // that was easy...
    dst[0] = scale[0] * pt[0] + offset[0];
    dst[1] = scale[1] * pt[1] + offset[1];
}

void draw_object(double (*verts)[3], size_t n_verts, RenderCtx *ctx)
{
    static double t = 0;

    vec3 centroid;
    mesh_centroid(centroid, verts, n_verts);
    vec3 camera_pos = {0, t+=0.001, 0};
    vec3 up = {0, 1, 0};
    mat4x4 camera_space_transform;
    camera_transform(camera_space_transform, camera_pos, verts[0], up);

    vec3 min_world_pt;
    vec3 max_world_pt;
    mesh_bounds(verts, n_verts, min_world_pt, max_world_pt);

    // removing leftover values from previous frame
    memset(ctx->buffer, 0, ctx->rows * ctx->cols * 3 * sizeof(uint8_t));
    memset(ctx->z_buffer, 0, ctx->rows * ctx->cols * sizeof(double));

    // looping through 3 verts at a time
    for (int i = 0; i < n_verts; i+=3)
    {
        vec3 screen_triangle[3];
        for (int j = 0; j < 3; j++)
        {
            vec4 norm_world_pt = {
                normalize(verts[i + j][0], max_world_pt[0], min_world_pt[0]),
                normalize(verts[i + j][1], max_world_pt[1], min_world_pt[1]),
                normalize(verts[i + j][2], max_world_pt[2], min_world_pt[2]),
                1
            };

            // 3D point in the camera's coordinate space
            vec4 camera_space_pt = {0};
            mat4x4_mul_vec4(camera_space_pt, camera_space_transform, norm_world_pt);

            vec2 clip_space_pt;
            vec2 scale = {1 * ctx->rows, 1 * ctx->cols};
            vec2 offset = {0,0};
            ortho_projection(clip_space_pt, camera_space_pt,scale,offset);

            // vec2, with the 3rd position used for the Z buffer
            vec3 screen_space_pt = {
                clip_space_pt[0],
                clip_space_pt[1],
                camera_space_pt[2]
            };

            memcpy(screen_triangle[j], screen_space_pt, sizeof(vec3));
        }

        vec3 normal;
        surface_normal(normal, &verts[i]);
        vec3_norm(normal, normal);

        uint8_t pixel[3] = {rand() % 255, rand() % 255, rand() % 255};
        draw_triangle(screen_triangle, pixel, ctx);
    }
}

RenderCtx init_renderer()
{
    RenderCtx ctx;
    ObjMesh *mesh = malloc(sizeof(ObjMesh));

    ctx.mesh = mesh;

    ctx.rows = 2048;
    ctx.cols = 2048;
    ctx.subpixel = 0.5;

    ctx.buffer = calloc(ctx.rows * ctx.cols * 3, sizeof(uint8_t));
    assert(ctx.buffer != NULL);

    ctx.z_buffer = calloc(ctx.rows * ctx.cols, sizeof(double));
    assert(ctx.z_buffer != NULL);

    // malloc'd by parse_obj
    parse_obj("models/teapot_maya_big.obj",
            &mesh->size, &mesh->verts, &mesh->texcoords, &mesh->normals);

    // parse_obj should return all faces as tris, but let's make sure
    assert(mesh->size % 3 == 0);

    return ctx;
}

void destroy_renderer(RenderCtx *ctx)
{
    free(ctx->mesh->verts);
    free(ctx->mesh->texcoords);
    free(ctx->mesh->normals);
    free(ctx->mesh);
    free(ctx->z_buffer);
    free(ctx->buffer);
}

