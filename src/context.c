#include "context.h"
#include "config.h"
#include "utils.h"
#include "obj_parser.h"
#include "rasterize.h"

#include <pthread.h>
#include <stdlib.h>
#include <assert.h>
#include <float.h>

static void divide_among_threads(int data_size, int n_threads, int out_sizes[n_threads]);
static void *clear_buffers_thread(void *v_ctx);

static void divide_among_threads(int data_size, int n_threads, int out_sizes[n_threads])
{
    memset(out_sizes, 0, sizeof(int[n_threads]));
    for (int i = 0; i < data_size; i++)
    {
        ++out_sizes[i % n_threads];
    }
}

static void *clear_buffers_thread(void *v_ctx)
{
    RenderCtx *ctx = v_ctx;

    memset(ctx->buffer_2, BACKGROUND, ctx->rows * ctx->cols * sizeof(color_t));

    // memset is much faster than a loop for some reason, 0x7f gets interpreted as a
    // high enough double value for the Z buffer test.
    memset(ctx->z_buffer_2, 0x7f, ctx->rows * ctx->cols * sizeof(double));

    return NULL;
}

// TODO: your threads are shit
pthread_t clear_buffers_start(RenderCtx *ctx)
{
    pthread_t thread = 0;
    clear_buffers_thread(ctx);
    /* int ret = pthread_create(&thread, NULL, clear_buffers_thread, ctx); */
    /* assert(ret == 0); */
    return thread;
}

void wait_for_clear(pthread_t thread)
{
    /* pthread_join(thread, NULL); */
}

void swap_buffers(RenderCtx *ctx)
{
    color_t (*tmp) = ctx->buffer;
    ctx->buffer = ctx->buffer_2;
    ctx->buffer_2 = tmp;

    double (*z_tmp) = ctx->z_buffer;
    ctx->z_buffer = ctx->z_buffer_2;
    ctx->z_buffer_2 = z_tmp;
}

void move_camera(RenderCtx *ctx, vec3 offset)
{
    mat4x4_translate_in_place(ctx->view_mat, -offset[0], -offset[1], -offset[2]);
}

RenderCtx init_renderer()
{
    RenderCtx ctx;
    Mesh *mesh = malloc(sizeof(Mesh));
    assert(mesh != NULL);

    ctx.mesh = mesh;

    ctx.rows = ROWS;
    ctx.cols = COLS;

    ctx.buffer = malloc(ctx.rows * ctx.cols * sizeof(color_t));
    ctx.buffer_2 = malloc(ctx.rows * ctx.cols * sizeof(color_t));
    assert(ctx.buffer != NULL && ctx.buffer_2 != NULL);

    ctx.z_buffer = malloc(ctx.rows * ctx.cols * sizeof(double));
    ctx.z_buffer_2 = malloc(ctx.rows * ctx.cols * sizeof(double));
    assert(ctx.z_buffer != NULL && ctx.z_buffer_2 != NULL);

    clear_buffers_thread(&ctx);

    parse_obj("./models/teapot_maya.obj",
            &mesh->size, &mesh->verts, &mesh->texcoords, &mesh->normals);
    // parse_obj should return all faces as tris, but let's make sure
    assert(mesh->size % 3 == 0);

    mesh_centroid(mesh->centroid, mesh->verts, mesh->size);

    ctx.grid_points = compute_grid(3, 5, &ctx.grid_rows, &ctx.grid_cols);

    vec3 up = UP_VECTOR;
    vec3 camera_pos = {0, 0, -2};
    /* vec3 cam_front = {0,0,1}; */
    /* vec3 target; */
    /* vec3_add(target, camera_pos, cam_front); */
    /* lookat(camera_pos, target, ctx.view_mat, up); */
    lookat(camera_pos, ctx.mesh->centroid, ctx.view_mat, up);

    // number of triangles per thread, need to multiply by 3 to get num points
    ctx.n_threads = 1;
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
    free(ctx->buffer);
    free(ctx->buffer_2);
    free(ctx->z_buffer);
    free(ctx->z_buffer_2);
    free(ctx->grid_points);
    free(ctx->thread_sizes);
}

