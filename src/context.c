#include "context.h"
#include "config.h"
#include "utils.h"
#include "obj_parser.h"
#include "rasterize.h"

#include <stdlib.h>
#include <assert.h>
#include <float.h>

void clear_buffers(RenderCtx *ctx)
{
    memset(ctx->buffer, BACKGROUND, ctx->rows * ctx->cols * sizeof(color_t));

    // memset is much faster than a loop for some reason, 0x7f gets interpreted as a
    // high enough double value for the Z buffer test.
    memset(ctx->z_buffer, 0x7f, ctx->rows * ctx->cols * sizeof(double));
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
    assert(ctx.buffer != NULL);

    ctx.z_buffer = malloc(ctx.rows * ctx.cols * sizeof(double));
    assert(ctx.z_buffer != NULL);

    clear_buffers(&ctx);

    parse_obj("./models/teapot_maya.obj",
            &mesh->size, &mesh->verts, &mesh->texcoords, &mesh->normals);
    // parse_obj should return all faces as tris, but let's make sure
    assert(mesh->size % 3 == 0);

    mesh_centroid(mesh->centroid, mesh->size, 3, mesh->verts);

    ctx.grid_points = compute_grid(3, 5, &ctx.grid_rows, &ctx.grid_cols);

    vec3 up = UP_VECTOR;
    vec3 camera_pos = {0, 0, -1};
    lookat(camera_pos, ctx.mesh->centroid, ctx.view_mat, up);

    return ctx;
}

void destroy_renderer(RenderCtx *ctx)
{
    free(ctx->mesh->verts);
    free(ctx->mesh->texcoords);
    free(ctx->mesh->normals);
    free(ctx->mesh);
    free(ctx->buffer);
    free(ctx->z_buffer);
    free(ctx->grid_points);
}

