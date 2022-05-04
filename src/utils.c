#include "utils.h"
#include "linmath_d.h"

#include <stdlib.h>
#include <assert.h>
#include <stdbool.h>
#include <float.h>

void print_vecn(double* vecn, int n)
{
    for (int i = 0; i < n; i++)
    {
        printf("%f, ", vecn[i]);
    }
    putchar('\n');
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

int ptr_dist(void *p1, void* p2)
{
    return (((char *)p2) - ((char *)p1));
}

int imin(int a, int b)
{
    return a > b ? b : a;
}

int imax(int a, int b)
{
    return a > b ? a : b;
}

bool dbl_almost_equal(double a, double b, double eps)
{
    return fabs(a - b) <= ((fabs(a) < fabs(b) ? fabs(b) : fabs(a)) * eps);
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

double lerp(double a, double b, double f)
{
    return (a * (1.0f - f)) + (b * f);
}

double inv_lerp(double a, double b, double x)
{
    return (x - a) / (b - a);
}

bool is_inbetween(double bound1, double bound2, double val)
{

    double lower = fmin(bound1, bound2);
    double upper = fmax(bound1, bound2);
    return (val < upper && val > lower);
}

bool point_in_bbox(vec2 p, double max_x, double min_x, double max_y, double min_y)
{
    return p[0] <= max_x && p[0] >= min_x && p[1] <= max_y && p[1] >= min_y;
}

// https://en.wikipedia.org/wiki/Line%E2%80%93plane_intersection
void plane_line_intersection(vec3 dst, vec3 a, vec3 b, vec3 plane_pt, vec3 plane_normal)
{
    vec3 tmp;
    vec3 direction;
    vec3_sub(direction, b, a);

    vec3_sub(tmp, plane_pt, a);
    double num = vec3_mul_inner(tmp, plane_normal);
    double denom = vec3_mul_inner(direction, plane_normal);

    // Not handling other cases, if this ever trips, i'll handle them.
    assert(num != 0 && denom != 0);

    double d = num/denom;

    vec3_scale(tmp, direction, d);
    vec3_add(dst, tmp, a);
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

void mesh_centroid(vec3 dst, int n_verts, int vert_size, double (*verts)[vert_size])
{
    assert(vert_size >= 3);
    long double sum_x, sum_y, sum_z;
    sum_x = sum_y = sum_z = 0;
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
void triangle_normal(vec3 dst, int tri_pt_size, double triangle[TRI_NPTS][tri_pt_size])
{
    vec3 v;
    vec3 w;
    vec3_sub(v, triangle[2], triangle[0]);
    vec3_sub(w, triangle[1], triangle[0]);
    vec3_mul_cross(dst, v, w);
    vec3_norm(dst, dst);
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

