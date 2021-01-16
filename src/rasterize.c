#include "rasterize.h"
#include "linmath_d.h"
#include "utils.h"
#include "config.h"
#include "context.h"
#include "draw.h"
#include "list.h"

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <float.h>
#include <stdbool.h>

enum ClipFlag{CLIPPED_TRI, CLIPPED_QUAD, UNCLIPPED, DONT_DRAW};

static void ortho_projection(vec2 dst, vec3 pt, vec2 scale, vec2 offset);
static void gl_projection_mat(mat4x4 ret, double t, double b, double l, double r, double n, double f);
static void view_frustum_bounds(double out_bounds[4], double angle_of_view, double aspect, double z);
static void world_to_clip(vec4 dst, vec4 world_pt, mat4x4 view_mat, double aspect, double fov_deg, double near, double far);
static void clip_space_to_screen(vec4 out_screen_pt, vec4 clip_space_pt, int rows, int cols);
static enum ClipFlag clip_triangle_near(vec4 triangle[TRI_NPTS], vec4 out_pts[TRI_NPTS * 2], double near_clip);
static bool cull_triangle(vec4 triangle[TRI_NPTS]);
static void rotate_point(vec3 dst, vec3 axis, vec3 position, double theta);
static double light_triangle(vec3 camera_z, vec3 world_tri[TRI_NPTS]);

void lookat(vec3 camera_pos, vec3 target, mat4x4 dst, vec3 up)
{
    // camera forward vector
    vec3 z_axis;
    vec3_sub(z_axis, camera_pos, target);
    vec3_norm(z_axis, z_axis);

    vec3 x_axis; // camera right vector
    vec3_mul_cross(x_axis, up, z_axis);

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
// TODO: This flips Y for some reason
static void ortho_projection(vec2 dst, vec3 pt, vec2 scale, vec2 offset)
{
    dst[0] = scale[0] * pt[0] + offset[0];
    dst[1] = scale[1] * pt[1] + offset[1];
}

// http://www.songho.ca/opengl/gl_projectionmatrix.html
static void gl_projection_mat(mat4x4 ret, double t, double b, double l, double r, double n, double f)
{
    ret[0][0] = 2 * n / (r - l);
    ret[0][1] = 0;
    ret[0][2] = 0;
    ret[0][3] = 0;

    ret[1][0] = 0;
    ret[1][1] = 2 * n / (t - b);
    ret[1][2] = 0;
    ret[1][3] = 0;

    ret[2][0] = (r + l) / (r - l);
    ret[2][1] = (t + b) / (t - b);
    ret[2][2] = -(f + n) / (f - n);
    ret[2][3] = -1;

    ret[3][0] = 0;
    ret[3][1] = 0;
    ret[3][2] = -2 * f * n / (f - n);
    ret[3][3] = 0;
}

// checks bounds of the view frustum at a given z value in view space
static void view_frustum_bounds(double out_bounds[4], double angle_of_view, double aspect, double z)
{
    double scale  = tan(angle_of_view * 0.5 * M_PI / 180) * z;
    double right  = aspect * scale;
    double left   = -right;
    double top    = scale;
    double bottom = -top;

    out_bounds[0] = top;
    out_bounds[1] = bottom;
    out_bounds[2] = left;
    out_bounds[3] = right;
}

static void world_to_clip(vec4 dst, vec4 world_pt, mat4x4 view_mat, double aspect, double fov_deg, double near, double far)
{
    double fbounds[4];
    mat4x4 proj_mat;
    vec4 camera_space_pt = {0,0,0,1};

    // getting frustum bounds at -near, which is the near clip plane in front of the camera.
    view_frustum_bounds(fbounds, fov_deg, aspect, -near);
    gl_projection_mat(proj_mat, fbounds[0], fbounds[1], fbounds[2], fbounds[3], near, far);
    mat4x4_mul_vec4(camera_space_pt, view_mat, world_pt);

    mat4x4_mul_vec4(dst, proj_mat, camera_space_pt);
}

static void clip_space_to_screen(vec4 out_screen_pt, vec4 clip_space_pt, int rows, int cols)
{
    assert(clip_space_pt[3] != 0);
    // perspective division by w normalizes pt to -1,1 NDC coords
    vec3 ndc_pt;
    ndc_pt[0] = clip_space_pt[0] / clip_space_pt[3];
    ndc_pt[1] = clip_space_pt[1] / clip_space_pt[3];
    ndc_pt[2] = clip_space_pt[2] / clip_space_pt[3];

    // scaling & offsetting origin to the middle of the framebuffer
    out_screen_pt[0] = ((ndc_pt[0] + 1) / 2 * cols);
    out_screen_pt[1] = ((1 - ndc_pt[1]) / 2 * rows);

    out_screen_pt[2] = ndc_pt[2];
    out_screen_pt[3] = 1 / clip_space_pt[3]; // depth inverse for persp correct interpolation
}

// clips the clip space triangle against the near clip plane
// clipped output is either one or two ordered triangles.
static enum ClipFlag clip_triangle_near(vec4 triangle[TRI_NPTS], vec4 out_pts[TRI_NPTS * 2], double near_clip)
{
    vec4 *clipped_pts[TRI_NPTS];
    vec4 *unclipped_pts[TRI_NPTS];
    int n_clipped_pts = 0;
    int n_unclipped_pts = 0;

    for (int i = 0; i < TRI_NPTS; i++)
    {
        if (triangle[i][3] < near_clip)
        {
            clipped_pts[n_clipped_pts++] = &triangle[i];
        }
        else
        {
            unclipped_pts[n_unclipped_pts++] = &triangle[i];
        }
    }

    // triangle doesn't touch clipping plane
    if (n_clipped_pts == 0)
    {
        return UNCLIPPED;
    }
    // triangle enveloped by clipping plane
    else if (n_unclipped_pts == 0)
    {
        return DONT_DRAW;
    }

    // 2 intersections between tri and clipping plane should be guaranteed from here

    // if there are 2 unclipped points, plus the 2 from the
    // plane intersections, the 4 points (quad) need to be triangulated.
    // otherwise, the return will only have 3 points (tri).
    bool clipped_quad = (n_unclipped_pts == 2);

    // When idx 1 is the only clipped point, the triangulated quad gets created with a clockwise (backwards) winding order
    // I'm not sure if this is always the case, but flipping it based on this seems to work.
    bool quad_flip = (clipped_quad) && ((int)(ptr_dist(triangle[0], clipped_pts[0]) / sizeof(vec4)) == 1);

    vec3 plane_pt = {0, 0, -near_clip}; // in front of camera, negative z
    vec3 plane_normal = {0, 0, 1};
    for (int i = 0; i < n_unclipped_pts; i++)
    {
        for (int j = 0; j < n_clipped_pts; j++)
        {
            vec4 intersect;

            intersect[3] = near_clip;
            plane_line_intersection(intersect, *unclipped_pts[i], *clipped_pts[j],
                    plane_pt, plane_normal);

            if (clipped_quad)
            {
                // first triangle  (i==0) -> {intersect, unclipped[i], unclipped[i + 1]}
                // second triangle (i==1) -> {unclipped[i], intersect, prev_intersect}

                if (i == 0)
                {
                    memcpy(out_pts[0], intersect, sizeof(vec4));
                    memcpy(out_pts[quad_flip ? 2 : 1], *unclipped_pts[i], sizeof(vec4));
                    memcpy(out_pts[quad_flip ? 1 : 2], *unclipped_pts[i + 1], sizeof(vec4));
                }
                else // i == 1
                {
                    memcpy(out_pts[3], *unclipped_pts[i], sizeof(vec4));
                    memcpy(out_pts[quad_flip ? 5 : 4], intersect, sizeof(vec4));
                    memcpy(out_pts[quad_flip ? 4 : 5], out_pts[0], sizeof(vec4));
                }

            }
            else // clipped_tri
            {
                // intersected point inherits it's clipped idx to keep winding order consistent
                int orig_idx = ptr_dist(triangle[0], clipped_pts[j]) / sizeof(vec4);
                memcpy(out_pts[orig_idx], intersect, sizeof(vec4));
            }
        }

        if (!clipped_quad)
        {
            int orig_idx = ptr_dist(triangle[0], unclipped_pts[i]) / sizeof(vec4);
            memcpy(out_pts[orig_idx], *unclipped_pts[i], sizeof(vec4));
        }
    }
    return clipped_quad ? CLIPPED_QUAD : CLIPPED_TRI;
}

// culls clip space triangle against viewport
static bool cull_triangle(vec4 triangle[TRI_NPTS])
{
    int clip_mask[TRI_NPTS];
    for (int i = 0; i < TRI_NPTS; i++)
    {
        double w = triangle[i][3];

        if (is_inbetween(w, -w, triangle[i][0]) &&
            is_inbetween(w, -w, triangle[i][1]))
        {
            clip_mask[i] = 0;
        }
        else
        {
            clip_mask[i] = 1;
        }
    }

    // only culling triangles fully outside viewport
    //      draw_triangle() takes care of clipping at screen edges
    return (clip_mask[0] && clip_mask[1] && clip_mask[2]);
}

// https://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula
static void rotate_point(vec3 dst, vec3 axis, vec3 position, double theta)
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

static double light_triangle(vec3 camera_z, vec3 world_tri[TRI_NPTS])
{
    vec3 normal;
    triangle_normal(normal, 3, world_tri);

    // light direction
    vec3 neg_cam = {-camera_z[0],-camera_z[1],-camera_z[2]};
    return vec3_mul_inner(normal, neg_cam);
}

void draw_object(RenderCtx *ctx)
{
    double time = inc(0.01, 360);

    const double aspect = (double) ctx->rows / ctx->cols;

    color_t tri_color = {90, 230, 95};

    for (int i = 0; i < ctx->mesh->size;  i+=TRI_NPTS)
    {

        vec4 clip_space_tri[TRI_NPTS];
        vec4 clipped_pts[TRI_NPTS * 2]; // clipped tri gets split into max 2 other tris

        for (int j = 0; j < TRI_NPTS; j++)
        {
            vec4 rotated_pt = {0, 0, 0, 1};
            vec3 norm_rot_axis = {0, 1, 0};
            rotate_point(rotated_pt, norm_rot_axis, ctx->mesh->verts[i + j], time);

            world_to_clip(clip_space_tri[j], rotated_pt, ctx->view_mat, aspect, FOV, NEAR_CLIP, FAR_CLIP);
        }

        // triangle is out of frustum bounds
        if (cull_triangle(clip_space_tri))
        {
            continue;
        }

        enum ClipFlag flag = clip_triangle_near(clip_space_tri, clipped_pts, NEAR_CLIP);

        switch (flag)
        {
            case UNCLIPPED:
            {
                vec4 screen_tri[TRI_NPTS];
                for (int i = 0; i < TRI_NPTS; i++)
                {
                    clip_space_to_screen(screen_tri[i], clip_space_tri[i], ctx->rows, ctx->cols);
                }
                draw_triangle(screen_tri, tri_color, false, ctx);
                /* for (int i = 0; i < TRI_NPTS; i++) */
                /* { */
                /*     draw_point(screen_tri[i], tri_color, 2, ctx); */
                /* } */
                break;
            }
            case CLIPPED_TRI:
            {
                vec4 screen_tri[TRI_NPTS];
                for (int i = 0; i < TRI_NPTS; i++)
                {
                    clip_space_to_screen(screen_tri[i], clipped_pts[i], ctx->rows, ctx->cols);
                }
                color_t t_color = {0, 0, 255};
                draw_triangle(screen_tri, t_color, false, ctx);
                break;
            }
            case CLIPPED_QUAD:
            {
                for (int i = 0; i < TRI_NPTS * 2; i+=TRI_NPTS)
                {
                    vec4 screen_tri[TRI_NPTS];
                    for (int j = 0; j < TRI_NPTS; j++)
                    {
                        clip_space_to_screen(screen_tri[j], clipped_pts[i + j], ctx->rows, ctx->cols);
                    }
                    color_t q_color = {255, 0, 0};
                    draw_triangle(screen_tri, q_color, false, ctx);
                }
                break;
            }
            case DONT_DRAW:
            {
                // Do nothing
                break;
            }
        }
    }
}

void draw_object_wireframe(RenderCtx *ctx)
{
    /* double time = inc(0.01, 1000); */
    /* mat4x4 view_mat; */
    /* lookat(ctx->camera_pos, ctx->mesh->centroid, */
    /*         view_mat, ctx->camera_z); */
    /* for (int i = 0; i < ctx->mesh->size; i += TRI_NPTS) */
    /* { */
    /*     vec3 screen_triangle[TRI_NPTS]; */
    /*     for (int j = 0; j < TRI_NPTS; j++) */
    /*     { */
    /*         vec4 rotated_pt={0, 0, 0, 1}; */
    /*         vec3 norm_rot_axis = {0, 1, 0}; */
    /*         rotate_point(rotated_pt, norm_rot_axis, ctx->mesh->verts[i + j], time); */

    /*         world_to_screen(screen_triangle[j], view_mat, rotated_pt, ctx->rows, ctx->cols); */
    /*     } */

    /*     color_t color = {102, 1, 19}; */
    /*     draw_line(screen_triangle[0][0], screen_triangle[0][1], screen_triangle[1][0], screen_triangle[1][1], */
    /*             color, ctx); */
    /*     draw_line(screen_triangle[1][0], screen_triangle[1][1], screen_triangle[2][0], screen_triangle[2][1], */
    /*             color, ctx); */
    /*     draw_line(screen_triangle[2][0], screen_triangle[2][1], screen_triangle[0][0], screen_triangle[0][1], */
    /*             color, ctx); */
    /* } */
}

vec2 **compute_grid(double size, int n_subdivs, int *out_rows, int *out_cols)
{
    // 4 lines in the main grid quad
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
    /* // storing grid_points as a double pointer, casting back to 2D array. */
    /* vec2 (*grid_pts_2d)[ctx->grid_cols] = (vec2 (*)[ctx->grid_cols])ctx->grid_points; */

    /* mat4x4 view_mat; */
    /* lookat(ctx->camera_pos, ctx->mesh->centroid, */
    /*         view_mat, ctx->camera_z); */
    /* for (int i = 0; i < ctx->grid_cols; i++) */
    /* { */
    /*     vec3 screen_pts[ctx->grid_rows]; */
    /*     for (int j = 0; j < ctx->grid_rows; j++) */
    /*     { */
    /*         vec4 world_pt = {0, 0, 0, 1}; */

    /*         world_pt[0] = grid_pts_2d[j][i][0]; */
    /*         world_pt[1] = 0; */
    /*         world_pt[2] = grid_pts_2d[j][i][1]; */

    /*         world_to_screen(screen_pts[j], view_mat, world_pt, ctx->rows, ctx->cols); */
    /*     } */

    /*     color_t color = {BACKGROUND-20, BACKGROUND-20, BACKGROUND-20}; */
    /*     if (i == ctx->grid_cols / 2) // Middle line */
    /*     { */
    /*         color[0] = 227; */
    /*         color[1] = 59; */
    /*         color[2] = 104; */
    /*     } */

    /*     draw_line(screen_pts[0][0], screen_pts[0][1], screen_pts[2][0], screen_pts[2][1], color, ctx); */
    /*     draw_line(screen_pts[1][0], screen_pts[1][1], screen_pts[3][0], screen_pts[3][1], color, ctx); */
    /* } */
}
