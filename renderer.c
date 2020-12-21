#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <stdint.h>
#include <assert.h>
#include <limits.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include <float.h>
#include <pthread.h>

// https://github.com/datenwolf/linmath.h
//   i made all the vec# types doubles, added string.h
#include "lib/linmath_d.h"

#include "obj_parser.h"
#include "wu_line.h"
#include "renderer.h"

#define ROWS 1024
#define COLS 1024
#define BACKGROUND 50
#define UP_VECTOR {0,1,0}
#define N_THREADS 48

void print_vec3(vec3 v)
{
    printf("%f %f %f \n", v[0], v[1], v[2]);
}

double inc(double i, double lim)
{
    static double t = 0;
    static bool sign = false;
    t += sign ? i : -i;
    if (t > lim)
        sign = false;
    else if (t < -lim)
        sign = true;
    return t;
}

int clamp(int val, int min, int max)
{
    if (val > max)
        return max;
    else if (val < min)
        return min;
    return val;
}

double normalize(double val, double upper, double lower)
{
    return (val - lower) / (upper - lower);
}

double gamma_correct(double val, double g)
{
    return 255. * pow((val / 255), g);
}

double distance(vec3 a, vec3 b)
{
    return sqrt(pow(b[0] - a[0], 2) + pow(b[1] - a[1], 2) + pow(b[2] - a[2], 2));
}

void line_midpoint(vec2 dst, vec2 a, vec2 b)
{
    double abs_mid_x = fabs(a[0] - b[0]) / 2.;
    double abs_mid_y = fabs(a[1] - b[1]) / 2.;
    dst[0] = fmin(a[0], b[0]) + abs_mid_x;
    dst[1] = fmin(a[1], b[1]) + abs_mid_y;
}

// Only fills the middle of the line, leaves the ends empty
void subdivide_line(vec2 *dst, size_t out_len, vec2 line[2], int n_subdivs)
{
    if (n_subdivs > 0)
    {
        size_t mid_idx = out_len / 2;

        vec2 mid;
        line_midpoint(mid, line[0], line[1]);
        memcpy(dst[mid_idx], mid, sizeof(vec2));

        vec2 subline_1[2];
        vec2 subline_2[2];
        memcpy(subline_1[0], line[0], sizeof(vec2));
        memcpy(subline_1[1], mid, sizeof(vec2));
        memcpy(subline_2[0], mid, sizeof(vec2));
        memcpy(subline_2[1], line[1], sizeof(vec2));

        size_t new_out_len = pow(2, n_subdivs - 1) + 1;
        subdivide_line(dst, new_out_len, subline_1, n_subdivs - 1);

        // assigning to the 2nd half of the output array
        subdivide_line(dst + mid_idx, new_out_len, subline_2, n_subdivs - 1);
    }
}

void mesh_bounds(vec3 *verts, size_t n_verts, vec3 out_min, vec3 out_max)
{
    double min_x, min_y, min_z;
    double max_x, max_y, max_z;
    min_x = min_y = min_z = DBL_MAX;
    max_x = max_y = max_z = -DBL_MAX;

    for (int i = 0; i < n_verts; i++)
    {
        if (verts[i][0] < min_x)
            min_x = verts[i][0];
        else if (verts[i][0] > max_x)
            max_x = verts[i][0];

        if (verts[i][1] < min_y)
            min_y = verts[i][1];
        else if (verts[i][1] > max_y)
            max_y = verts[i][1];

        if (verts[i][2] < min_z)
            min_z = verts[i][2];
        else if (verts[i][2] > max_z)
            max_z = verts[i][2];
    }

    // everything should have a value
    assert(min_x != DBL_MAX && min_y != DBL_MAX && min_z != DBL_MAX);
    assert(max_x != -DBL_MAX && max_y != -DBL_MAX && max_z != -DBL_MAX);

    vec3 mins = {min_x, min_y, min_z};
    vec3 maxs = {max_x, max_y, max_z};
    memcpy(out_min, mins, sizeof(vec3));
    memcpy(out_max, maxs, sizeof(vec3));
}

void mesh_centroid(vec3 dst, vec3 *verts, size_t n_verts)
{
    long double sum_x = 0, sum_y = 0, sum_z = 0;
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

// computes surface normal of triangle
void surface_normal(vec3 dst, vec3 *triangle)
{
    vec3 v;
    vec3 w;
    vec3_sub(v, triangle[2], triangle[0]);
    vec3_sub(w, triangle[1], triangle[0]);
    vec3_mul_cross(dst, v, w);
}

void divide_among_threads(int data_size, int n_threads, int out_sizes[n_threads])
{
    memset(out_sizes, 0, sizeof(int[n_threads]));
    for (int i = 0; i < data_size; i++)
    {
        ++out_sizes[i % n_threads];
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

void normalize_coords(vec3 *out_arr, vec3 *verts, size_t n_verts)
{
    vec3 min_pt;
    vec3 max_pt;
    mesh_bounds(verts, n_verts, min_pt, max_pt);

    // I'm not sure if this is right, but I think normalizing
    // each coord to it's respective max/min value would skew the mesh.
    // So i'm normalizing all of the coords to the same range
    double max_max_pt = fmax(max_pt[0],fmax(max_pt[1],max_pt[2]));
    double min_min_pt = fmin(min_pt[0],fmin(min_pt[1],min_pt[2]));

    for (int i = 0; i < n_verts; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            // normalizing between -1 and 1
            out_arr[i][j] = (normalize(verts[i][j], max_max_pt, min_min_pt) - 0.5) * 2.;
        }
    }
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
void draw_triangle(vec3 triangle[3], uint8_t fill[3], RenderCtx *ctx)
{
    double *a = triangle[0];
    double *b = triangle[1];
    double *c = triangle[2];

    double max_x, max_y, min_x, min_y;
    triangle_bbox(triangle, &max_x, &max_y, &min_x, &min_y);

    max_x = ceil(max_x);
    max_y = ceil(max_y);
    min_y = floor(min_y);
    min_x = floor(min_x);

    vec3 sp; // screen space point
    for (sp[0] = min_x; sp[0] <= max_x; sp[0] += 1)
    {
        for (sp[1] = min_y; sp[1] <= max_y; sp[1] += 1)
        {
            // if bbox point is out of screen bounds
            if (sp[0] >= ctx->cols - 1 || sp[1] >= ctx->rows - 1 || sp[0] < 0 || sp[1] < 0)
            {
                continue;
            }

            // barycentric coords
            double w0 = orient2d(b, c, sp);
            double w1 = orient2d(c, a, sp);
            double w2 = orient2d(a, b, sp);

            // if point is inside or on all edges
            if (w0 >= 0 && w1 >= 0 && w2 >= 0)
            {
                // init value for Z buffer
                sp[2] = 0;

                // smearing actual z values over bayercentric coords
                // (This only works for orthographic projections)
                sp[2] += a[2] * w0;
                sp[2] += b[2] * w1;
                sp[2] += c[2] * w2;

                // flipping vertically
                int row = (ctx->rows - 1) - sp[1];
                int col = sp[0];

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
    memset(ctx->buffer, BACKGROUND, ctx->rows * ctx->cols * 3 * sizeof(uint8_t));

    // This is much faster than a loop, -10000 is (hopefully) casted to a large negative double
    /* memset(ctx->z_buffer, -10000, ctx->rows * ctx->cols  * sizeof(double)); */

    // TODO: this is really slow
    for (int i = 0; i < ctx->rows * ctx->cols; ++i)
    {
        ctx->z_buffer[i] = -DBL_MAX;
    }
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
        uint8_t c = gamma_correct(clamp(intensity * 255, 0, 255), 2.2);
        uint8_t tri_color[3] = {c, c, c};

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

        uint8_t color[3] = {180, 255, 255};
        draw_line(screen_triangle[0][0], screen_triangle[0][1], screen_triangle[1][0], screen_triangle[1][1],
                ctx->rows, ctx->cols, ctx->buffer, color);
        draw_line(screen_triangle[1][0], screen_triangle[1][1], screen_triangle[2][0], screen_triangle[2][1],
                ctx->rows, ctx->cols, ctx->buffer, color);
        draw_line(screen_triangle[2][0], screen_triangle[2][1], screen_triangle[0][0], screen_triangle[0][1],
                ctx->rows, ctx->cols, ctx->buffer, color);
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

        uint8_t color[3] = {BACKGROUND-20, BACKGROUND-20, BACKGROUND-20};
        if (i == ctx->grid_cols / 2) // Middle line
        {
            color[0] = 227;
            color[1] = 59;
            color[2] = 104;
        }

        draw_line(screen_pts[0][0], screen_pts[0][1], screen_pts[2][0], screen_pts[2][1],
                ctx->rows, ctx->cols, ctx->buffer, color);
        draw_line(screen_pts[1][0], screen_pts[1][1], screen_pts[3][0], screen_pts[3][1],
                ctx->rows, ctx->cols, ctx->buffer, color);
    }
}

RenderCtx init_renderer()
{
    RenderCtx ctx;
    ObjMesh *mesh = malloc(sizeof(ObjMesh));
    assert(mesh != NULL);

    ctx.mesh = mesh;

    ctx.rows = ROWS;
    ctx.cols = COLS;

    ctx.buffer = calloc(ctx.rows * ctx.cols * 3, sizeof(uint8_t));
    assert(ctx.buffer != NULL);

    ctx.z_buffer = calloc(ctx.rows * ctx.cols, sizeof(double));
    assert(ctx.z_buffer != NULL);

    parse_obj("./models/teapot_maya.obj",
            &mesh->size, &mesh->verts, &mesh->texcoords, &mesh->normals);
    // parse_obj should return all faces as tris, but let's make sure
    assert(mesh->size % 3 == 0);

    // normalizing mesh coordinates in place
    normalize_coords(mesh->verts, mesh->verts, mesh->size);

    mesh_centroid(mesh->centroid, mesh->verts, mesh->size);

    ctx.camera_pos[0] = 0;
    ctx.camera_pos[1] = 0;
    ctx.camera_pos[2] = 0;

    vec3 up = UP_VECTOR;
    camera_transform(ctx.camera_pos, ctx.mesh->centroid, up,
            ctx.cam_transform, ctx.camera_z);

    ctx.grid_points = compute_grid(3, 5, &ctx.grid_rows, &ctx.grid_cols);

    // number of triangles per thread, need to multiply by 3 to get num points
    ctx.n_threads = N_THREADS;
    ctx.thread_sizes = malloc(sizeof(int) * ctx.n_threads);
    divide_among_threads(ctx.mesh->size / 3, ctx.n_threads, ctx.thread_sizes);

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
    free(ctx->grid_points);
    free(ctx->thread_sizes);
}

