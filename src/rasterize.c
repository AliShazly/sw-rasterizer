#include "rasterize.h"
#include "linmath_d.h"
#include "utils.h"
#include "config.h"
#include "context.h"
#include "draw.h"

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <float.h>
#include <pthread.h>

typedef struct
{
    size_t start_idx;
    size_t end_idx;
    double time; // temporary
    RenderCtx *ctx;

}ThreadArgs;

// look-at
void camera_transform(vec3 camera_pos, vec3 target, vec3 up, mat4x4 dst, vec3 out_z)
{
    // keeping target point at y=0
    /* vec3 t = {target[0], 0, target[1]}; */

    // camera forward vector
    vec3_sub(out_z, target, camera_pos);
    vec3_norm(out_z, out_z);

    vec3 x_axis; // camera right vector
    vec3_mul_cross(x_axis, up, out_z);
    vec3_norm(x_axis, x_axis);

    vec3 y_axis; // camera up vector
    vec3_mul_cross(y_axis, out_z, x_axis);

    /* mat4x4 orientation; */
    /* mat4x4_identity(orientation); */
    mat4x4 orientation = {
        {x_axis[0], y_axis[0], out_z[0], 0},
        {x_axis[1], y_axis[1], out_z[1], 0},
        {x_axis[2], y_axis[2], out_z[2], 0},
        {0,         0,         0,        1}
    };

    mat4x4 translation;
    mat4x4_identity(translation);
    translation[3][0] = -camera_pos[0];
    translation[3][1] = -camera_pos[1];
    translation[3][2] = -camera_pos[2];

    mat4x4_mul(dst, orientation, translation);
}

// https://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula
void rotate_point(vec3 dst, vec3 axis, vec3 position, double theta)
{
    vec3 d;
    vec3_scale(d, axis, vec3_mul_inner(axis, position));

    vec3 r;
    vec3_sub(r, position, d);

    vec3 rp1;
    vec3 rp2;
    vec3 rp3;
    vec3 rp;
    vec3_scale(rp1, r, cos(theta));
    vec3_mul_cross(rp2, axis, position);
    vec3_scale(rp3, rp2, sin(theta));
    vec3_add(rp, rp1, rp3);

    vec3_add(dst, d, rp);
}

// https://en.wikipedia.org/wiki/3D_projection#Orthographic_projection
void ortho_projection(vec2 dst, vec3 pt, vec2 scale, vec2 offset)
{
    dst[0] = scale[0] * pt[0] + offset[0];
    dst[1] = scale[1] * pt[1] + offset[1];
}

void world_to_screen(vec3 dst, mat4x4 cam_space_transform, vec4 world_pt, int rows, int cols)
{
    // 3D point in the camera's coordinate space
    vec4 camera_space_pt = {0};
    mat4x4_mul_vec4(camera_space_pt, cam_space_transform, world_pt);

    // Camera point projected & scaled to screen coordinates
    vec2 clip_space_pt;
    vec2 scale = {0.3 * cols, 0.3 * rows};
    vec2 offset = {cols/2., rows/2.};
    ortho_projection(clip_space_pt, camera_space_pt, scale, offset);

    // vec2, with the 3rd position used for the Z buffer
    vec3 screen_space_pt = {
        clip_space_pt[0],
        clip_space_pt[1],
        camera_space_pt[2]
    };

    memcpy(dst, screen_space_pt, sizeof(vec3));
}

void clear_buffers(RenderCtx *ctx)
{
    memset(ctx->buffer, BACKGROUND, ctx->rows * ctx->cols * sizeof(color_t));

    // This is much faster than a loop, -10000 is (hopefully) casted to a large negative double
    /* memset(ctx->z_buffer, -10000, ctx->rows * ctx->cols  * sizeof(double)); */

    // TODO: this is really slow
    for (int i = 0; i < ctx->rows * ctx->cols; ++i)
    {
        ctx->z_buffer[i] = -DBL_MAX;
    }
}

double light_triangle(vec3 camera_z, vec3 world_tri[3])
{
    vec3 normal;
    surface_normal(normal, world_tri);
    vec3_norm(normal, normal);

    // light direction
    vec3 neg_cam = {-camera_z[0],-camera_z[1],-camera_z[2]};
    return vec3_mul_inner(normal, neg_cam);
}

void *object_thread(void *thread_args)
{
    ThreadArgs *args = thread_args;

    for (int i = args->start_idx; i < args->end_idx + 1;  i+=3)
    {
        vec3 screen_tri[3];
        vec3 world_tri[3];
        for (int j = 0; j < 3; j++)
        {
            vec4 rotated_pt={0, 0, 0, 1};
            vec3 norm_rot_axis = {0, 1, 0};
            rotate_point(rotated_pt, norm_rot_axis, args->ctx->mesh->verts[i + j], args->time);

            world_to_screen(screen_tri[j], args->ctx->cam_transform,
                    rotated_pt, args->ctx->rows, args->ctx->cols);
            memcpy(world_tri[j], rotated_pt, sizeof(vec3));
        }

        double intensity = light_triangle(args->ctx->camera_z, world_tri);
        int c = gamma_correct(clamp(intensity * 255, 0, 255), 2.2);
        color_t tri_color = {c, c, c};

        // not drawing faces looking away from camera (backface culling)
        if (intensity >= 0)
        {
            draw_triangle(screen_tri, tri_color, args->ctx);
        }
    }
    return NULL;
}

void draw_object_threads(RenderCtx *ctx)
{
    double time = inc(0.01, 1000);
    pthread_t threads[ctx->n_threads];
    ThreadArgs *args[ctx->n_threads];

    int offset = -1;
    // drawing each thread's triangles
    // TODO: this assumes all threads have an assigned size > 0
    for (int i = 0; i < ctx->n_threads; i++)
    {
        args[i] = malloc(sizeof(ThreadArgs));
        assert(args[i] != NULL);

        args[i]->start_idx = offset + 1;
        args[i]->end_idx = offset + (3 * ctx->thread_sizes[i]);
        offset = args[i]->end_idx;

        args[i]->time = time;

        args[i]->ctx = ctx;

        int ret = pthread_create(&threads[i], NULL, object_thread, args[i]);
        assert(ret == 0);
    }

    // the end_idx of the last thread should be the end of the vert array
    assert(offset == ctx->mesh->size - 1);

    for (int i = 0; i < ctx->n_threads; i++)
    {
        pthread_join(threads[i], NULL);
        free(args[i]);
    }
}

void draw_object_wireframe(RenderCtx *ctx)
{
    double time = inc(0.01, 1000);
    for (int i = 0; i < ctx->mesh->size; i += 3)
    {
        vec3 screen_triangle[3];
        for (int j = 0; j < 3; j++)
        {
            vec4 rotated_pt={0, 0, 0, 1};
            vec3 norm_rot_axis = {0, 1, 0};
            rotate_point(rotated_pt, norm_rot_axis, ctx->mesh->verts[i + j], time);

            world_to_screen(screen_triangle[j], ctx->cam_transform, rotated_pt, ctx->rows, ctx->cols);
        }

        color_t color = {180, 255, 255};
        draw_line(screen_triangle[0][0], screen_triangle[0][1], screen_triangle[1][0], screen_triangle[1][1],
                color, ctx);
        draw_line(screen_triangle[1][0], screen_triangle[1][1], screen_triangle[2][0], screen_triangle[2][1],
                color, ctx);
        draw_line(screen_triangle[2][0], screen_triangle[2][1], screen_triangle[0][0], screen_triangle[0][1],
                color, ctx);
    }
}

vec2 **compute_grid(double size, int n_subdivs, int *out_rows, int *out_cols)
{
    const int n_lines = 4;

    vec2 lines[n_lines][2] = {
        {{-size,  size}, {size,   size}}, // set 1
        {{size,   size}, {size,  -size}}, // set 2
        {{-size, -size}, {size,  -size}}, // set 1
        {{-size,  size}, {-size, -size}}  // set 2
    };

    size_t subd_line_len = pow(2, n_subdivs) + 1;
    vec2 (*out_subd_lines)[subd_line_len] = malloc(sizeof(vec2) * subd_line_len * n_lines);
    assert(out_subd_lines != NULL);

    for (int i = 0; i < n_lines; i++)
    {
        // filling the ends of the line array with the parent line values
        memcpy(out_subd_lines[i][0                ], lines[i][0], sizeof(vec2));
        memcpy(out_subd_lines[i][subd_line_len - 1], lines[i][1], sizeof(vec2));

        subdivide_line(out_subd_lines[i], subd_line_len, lines[i], n_subdivs);
    }

    *out_rows = n_lines;
    *out_cols = subd_line_len;
    return (vec2**)out_subd_lines; // need to cast back to 2D array
}

void draw_grid(RenderCtx *ctx)
{
    // storing grid_points as a double pointer, casting back to 2D array.
    vec2 (*grid_pts_2d)[ctx->grid_cols] = (vec2 (*)[ctx->grid_cols])ctx->grid_points;

    for (int i = 0; i < ctx->grid_cols; i++)
    {
        vec3 screen_pts[ctx->grid_rows];
        for (int j = 0; j < ctx->grid_rows; j++)
        {
            vec4 world_pt = {0, 0, 0, 1};

            world_pt[0] = grid_pts_2d[j][i][0];
            world_pt[1] = 0;
            world_pt[2] = grid_pts_2d[j][i][1];
            world_to_screen(screen_pts[j], ctx->cam_transform, world_pt, ctx->rows, ctx->cols);
        }

        color_t color = {BACKGROUND-20, BACKGROUND-20, BACKGROUND-20};
        if (i == ctx->grid_cols / 2) // Middle line
        {
            color[0] = 227;
            color[1] = 59;
            color[2] = 104;
        }

        draw_line(screen_pts[0][0], screen_pts[0][1], screen_pts[2][0], screen_pts[2][1], color, ctx);
        draw_line(screen_pts[1][0], screen_pts[1][1], screen_pts[3][0], screen_pts[3][1], color, ctx);
    }
}
