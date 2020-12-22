#include "context.h"
#include "config.h"
#include "utils.h"
#include "obj_parser.h"
#include "rasterize.h"

#include <stdlib.h>
#include <assert.h>

void divide_among_threads(int data_size, int n_threads, int out_sizes[n_threads])
{
    memset(out_sizes, 0, sizeof(int[n_threads]));
    for (int i = 0; i < data_size; i++)
    {
        ++out_sizes[i % n_threads];
    }
}

RenderCtx init_renderer()
{
    RenderCtx ctx;
    Mesh *mesh = malloc(sizeof(Mesh));
    assert(mesh != NULL);

    ctx.mesh = mesh;

    ctx.rows = ROWS;
    ctx.cols = COLS;

    ctx.buffer = calloc(ctx.rows * ctx.cols, sizeof(color_t));
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
    // TODO: if the camera moves, the transform needs to be recomputed
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

