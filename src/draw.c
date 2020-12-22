#include "draw.h"
#include "context.h"
#include "linmath_d.h"
#include "utils.h"

#include <stdlib.h>
#include <math.h>

// barycentric coords
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
void draw_triangle(vec3 triangle[3], color_t fill, RenderCtx *ctx)
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

    int sp[3]; // screen space point
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
            vec2 sp_v = {sp[0], sp[1]};
            double w0 = orient2d(b, c, sp_v);
            double w1 = orient2d(c, a, sp_v);
            double w2 = orient2d(a, b, sp_v);

            // if point is inside or on all edges
            if (w0 >= 0 && w1 >= 0 && w2 >= 0)
            {
                sp[2] = 0; // init value for z buffer

                // smearing actual z values over bayercentric coords
                // (This only works for orthographic projections)
                sp[2] += a[2] * w0;
                sp[2] += b[2] * w1;
                sp[2] += c[2] * w2;

                int idx = coord_to_idx(sp[0], sp[1], ctx->rows, ctx->cols);

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

// https://en.wikipedia.org/wiki/Bresenham%27s_line_algorithm#All_cases
void draw_line(int x0, int y0, int x1, int y1, color_t color, RenderCtx *ctx)
{
    int dx = abs(x1 - x0);
    int sx = x0 < x1 ? 1 : -1;
    int dy = -abs(y1 - y0);
    int sy = y0 < y1 ? : -1;
    int err = dx + dy;
    for(;;)
    {
        // if point is in bounds
        if ((x0 < ctx->cols && x0 > 0) && (y0 < ctx->rows && y0 > 0))
        {
            // plotting point
            int idx = coord_to_idx(x0, y0, ctx->rows, ctx->cols);
            for (int i = 0; i < 3; i++)
            {
                ctx->buffer[idx][i] = color[i];
            }
        }

        if (x0 == x1 && y0 == y1)
            break;
        int e2 = 2 * err;
        if (e2 > dy)
        {
            err += dy;
            x0 += sx;
        }
        if (e2 <= dx)
        {
            err += dx;
            y0 += sy;
        }
    }
}

