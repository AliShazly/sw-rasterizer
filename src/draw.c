#include "draw.h"
#include "fixed_point_32.h"
#include "linmath_d.h"
#include "context.h"
#include "utils.h"

#include <assert.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <float.h>
#include <stdint.h>

typedef fixed_t vec2_fixed[2];

static int coord_to_idx(int x, int y, int rows, int cols);
static bool top_left_edge_cclock(vec2_fixed p1, vec2_fixed p2);
static fixed_t edge_func(vec2_fixed a, vec2_fixed b, vec2_fixed c);
static bool point_in_circle(vec2 pt, vec2 center, int radius);
static void circle_bbox(vec2 center, int radius, int *max_x, int *max_y, int *min_x, int *min_y);


// x,y screen coord to corresponding 1D buffer array index
static int coord_to_idx(int x, int y, int rows, int cols)
{
    // flipping vertically
    int row = (rows - 1) - y;
    // 2D idx (row, col) to a 1D index
    return row * cols + x;
}

static bool top_left_edge_cclock(vec2_fixed p1, vec2_fixed p2)
{
    // left edge goes down, top edge is exactly horizontal and goes towards the left
    return ((p1[1] > p2[1]) || (p1[1] == p2[1] && p1[0] < p2[0]));
}

// point c in relation to edge ab
static fixed_t edge_func(vec2_fixed a, vec2_fixed b, vec2_fixed c)
{
    fixed_t m0 = fixed_mult((c[0] - a[0]), (b[1] - a[1]));
    fixed_t m1 = fixed_mult((c[1] - a[1]), (b[0] - a[0]));
    return m0 - m1;
}

static bool point_in_circle(vec2 pt, vec2 center, int radius)
{
    // no need for a square root, square the radius for comparison instead
    return (pow((center[0] - pt[0]), 2) + pow((center[1] - pt[1]), 2)) < pow(radius, 2);
}

static void circle_bbox(vec2 center, int radius, int *max_x, int *max_y, int *min_x, int *min_y)
{
    *max_x = center[0] + radius;
    *max_y = center[1] + radius;
    *min_x = center[0] - radius;
    *min_y = center[1] - radius;
}

// http://people.csail.mit.edu/ericchan/bib/pdf/p17-pineda.pdf
// this assumes triangle has a counter clockwise winding order
void draw_triangle(vec4 triangle[TRI_NPTS], color_t fill, bool two_sided, RenderCtx *ctx)
{
    // converting x/y coords to 26.6 fixed point
    vec2_fixed a = {
        fixed_from_dbl(triangle[0][0]),
        fixed_from_dbl(triangle[0][1])};
    vec2_fixed b = {
        fixed_from_dbl(triangle[1][0]),
        fixed_from_dbl(triangle[1][1])};
    vec2_fixed c = {
        fixed_from_dbl(triangle[2][0]),
        fixed_from_dbl(triangle[2][1])};

    fixed_t area = abs(edge_func(a, b, c));
    if (area == 0) // triangle is degenerate
    {
        return;
    }

    // reassigning the z and inv_w coords
    double a_z, b_z, c_z, a_inv_w, b_inv_w, c_inv_w;
    a_z = triangle[0][2], a_inv_w = triangle[0][3];
    b_z = triangle[1][2], b_inv_w = triangle[1][3];
    c_z = triangle[2][2], c_inv_w = triangle[2][3];

    // triangle bounding box
    int max_x, max_y, min_x, min_y;
    max_x = fixed_to_dbl(imax(a[0], imax(b[0], c[0])));
    max_y = fixed_to_dbl(imax(a[1], imax(b[1], c[1])));
    min_x = fixed_to_dbl(imin(a[0], imin(b[0], c[0])));
    min_y = fixed_to_dbl(imin(a[1], imin(b[1], c[1])));

    // clipping bbox to screen bounds
    max_x = imin(max_x, ctx->cols - 1);
    max_y = imin(max_y, ctx->rows - 1);
    min_y = imax(min_y, 0);
    min_x = imax(min_x, 0);

    // checking fill rules for each edge
    int bias_w0 = top_left_edge_cclock(b, c) ? 0 : 1;
    int bias_w1 = top_left_edge_cclock(c, a) ? 0 : 1;
    int bias_w2 = top_left_edge_cclock(a, b) ? 0 : 1;

    vec2_fixed min_corner = {
        fixed_from_dbl(min_x + 0.5),
        fixed_from_dbl(min_y + 0.5)};
    fixed_t w0_row = edge_func(b, c, min_corner) + fixed_from_dbl(bias_w0);
    fixed_t w1_row = edge_func(c, a, min_corner) + fixed_from_dbl(bias_w1);
    fixed_t w2_row = edge_func(a, b, min_corner) + fixed_from_dbl(bias_w2);

    const int pstep = 1; // pixel step
    // used for incrementing edge function output
    fixed_t w0_step_x = pstep * (c[1] - b[1]);
    fixed_t w0_step_y = pstep * (b[0] - c[0]);
    fixed_t w1_step_x = pstep * (a[1] - c[1]);
    fixed_t w1_step_y = pstep * (c[0] - a[0]);
    fixed_t w2_step_x = pstep * (b[1] - a[1]);
    fixed_t w2_step_y = pstep * (a[0] - b[0]);

    for (int sp_y = min_y; sp_y <= max_y; sp_y+=pstep)
    {
        fixed_t w0 = w0_row;
        fixed_t w1 = w1_row;
        fixed_t w2 = w2_row;

        for (int sp_x = min_x; sp_x <= max_x; sp_x+=pstep)
        {
            bool backface = (w0 < 0 && w1 < 0 && w2 < 0) && two_sided;
            if ((w0 > 0 && w1 > 0 && w2 > 0) || backface)
            {
                double w0_a = fixed_to_dbl(fixed_divide(w0, area));
                double w1_a = fixed_to_dbl(fixed_divide(w1, area));
                double w2_a = fixed_to_dbl(fixed_divide(w2, area));

                // interpolating z values over barycentric coords
                double pixel_z = (a_z * w0_a) + (b_z * w1_a) + (c_z * w2_a);
                double denom = (a_inv_w * w0_a) + (b_inv_w * w1_a) + (c_inv_w * w2_a);
                pixel_z /= denom;

                int buf_idx = coord_to_idx(sp_x, sp_y, ctx->rows, ctx->cols);

                // depth buffer test & drawing pixels
                if (pixel_z < ctx->z_buffer[buf_idx])
                {
                    ctx->z_buffer[buf_idx] = pixel_z;

                    ctx->buffer[buf_idx][0] = fill[0];
                    /* ctx->buffer[buf_idx][1] = fill[1]; */
                    ctx->buffer[buf_idx][2] = fill[2];

                    /* ctx->buffer[buf_idx][0] = w0_a * 255; */
                    ctx->buffer[buf_idx][1] = w1_a * 255;
                    /* ctx->buffer[buf_idx][2] = w2_a * 255; */

                    int p_v = clamp(normalize(pixel_z, 2.5, -.1) * 255, 0, 255);
                    /* ctx->buffer[buf_idx][0] = p_v; */
                    /* ctx->buffer[buf_idx][1] = p_v; */
                    /* ctx->buffer[buf_idx][2] = p_v; */
                }
            }
            w0 += w0_step_x;
            w1 += w1_step_x;
            w2 += w2_step_x;
        }
        w0_row += w0_step_y;
        w1_row += w1_step_y;
        w2_row += w2_step_y;
    }
}


void draw_point(vec2 p, color_t color, int radius, RenderCtx *ctx)
{
    int max_x, max_y, min_x, min_y;
    circle_bbox(p, radius, &max_x, &max_y, &min_x, &min_y);

    // clipping bbox to screen bounds
    max_x = imin(max_x, ctx->cols - 1);
    max_y = imin(max_y, ctx->rows - 1);
    min_y = imax(min_y, 0);
    min_x = imax(min_x, 0);

    vec2 sp;
    for (sp[0] = min_x; sp[0] <= max_x; sp[0]++)
    {
        for (sp[1] = min_y; sp[1] <= max_y; sp[1]++)
        {
            if (point_in_circle(sp, p, radius))
            {
                int idx = coord_to_idx(sp[0], sp[1], ctx->rows, ctx->cols);

                memcpy(ctx->buffer[idx], color, sizeof(color_t));
            }
        }
    }
}

void draw_line(vec2 p1, vec2 p2, color_t color, int thick, RenderCtx *ctx)
{
    vec2 thickness = {thick, thick};

    vec2 a2;
    vec2_add(a2, p1, thickness);

    vec2 b2;
    vec2_add(b2, p2, thickness);

    vec3 tri_a[3] = {
        {p1[0], p1[1], 1},
        {a2[0], a2[1], 1},
        {p2[0], p2[1], 1}
    };

    vec3 tri_b[3] = {
        {p2[0], p2[1], 1},
        {b2[0], b2[1], 1},
        {a2[0], a2[0], 1}
    };

    /* draw_triangle(tri_a, color, ctx); */
    /* draw_triangle(tri_b, color, ctx); */

}

// https://en.wikipedia.org/wiki/Bresenham%27s_line_algorithm#All_cases
/* void draw_line(int x0, int y0, int x1, int y1, color_t color, RenderCtx *ctx) */
/* { */
/*     int dx = abs(x1 - x0); */
/*     int sx = x0 < x1 ? 1 : -1; */
/*     int dy = -abs(y1 - y0); */
/*     int sy = y0 < y1 ? : -1; */
/*     int err = dx + dy; */
/*     for(;;) */
/*     { */
/*         // if point is in bounds */
/*         if ((x0 < ctx->cols && x0 > 0) && (y0 < ctx->rows && y0 > 0)) */
/*         { */
/*             // plotting point */
/*             int idx = coord_to_idx(x0, y0, ctx->rows, ctx->cols); */
/*             for (int i = 0; i < 3; i++) */
/*             { */
/*                 ctx->buffer[idx][i] = color[i]; */
/*             } */
/*         } */

/*         if (x0 == x1 && y0 == y1) */
/*             break; */
/*         int e2 = 2 * err; */
/*         if (e2 > dy) */
/*         { */
/*             err += dy; */
/*             x0 += sx; */
/*         } */
/*         if (e2 <= dx) */
/*         { */
/*             err += dx; */
/*             y0 += sy; */
/*         } */
/*     } */
/* } */

