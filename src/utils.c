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
    assert(upper > lower);
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

// https://www.mathworks.com/matlabcentral/answers/16243-angle-between-two-vectors-in-3d
double angle3D(vec3 a, vec3 b)
{
    vec3 c;
    vec3_mul_cross(c, a, b);
    double d = vec3_mul_inner(a, b);
    return atan2(vec3_magnitude(c), d);
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

    // Not handling other cases
    assert(num != 0 && denom != 0);

    double d = num/denom;

    vec3 intersection_pt;
    vec3_scale(tmp, direction, d);
    vec3_add(intersection_pt, tmp, a);
    memcpy(dst, intersection_pt, sizeof(vec3));
}

int orient_compare(const void *vc1, const void* vc2)
{
    // only need the first element of each struct for comparison
    double theta1 = *((double *)vc1);
    double theta2 = *((double *)vc2);

    if (theta1 > theta2)
        return 1;
    else if (theta1 < theta2)
        return -1;
    else
        // not worried about strict fp equality
        return 0;
}

void orient_convex_polygon(vec3 *poly_pts, int poly_len, bool clockwise)
{
    typedef struct
    {
        double angle_to_centroid;
        int    pt_index;
    }VertCompare;

    vec3 centroid;
    mesh_centroid(centroid, poly_pts, poly_len);

    VertCompare poly_compare[poly_len];
    const int scale = clockwise ? 1 : -1; // sorting backwards for c-clockwise
    for (int i = 0; i < poly_len; i++)
    {
        poly_compare[i].angle_to_centroid = angle3D(poly_pts[i], centroid) * scale;
        poly_compare[i].pt_index = i;
    }

    // this fails if all angles are the same
    qsort(poly_compare, poly_len, sizeof(VertCompare), orient_compare);

    vec3 tmp_pts[poly_len];
    memcpy(tmp_pts, poly_pts, poly_len * sizeof(vec3));
    for (int i = 0; i < poly_len; i++)
    {
        int new_idx = poly_compare[i].pt_index;
        memcpy(poly_pts[new_idx], tmp_pts[i], sizeof(vec3));
    }
}

int fan_triangulate_out_len(int poly_len)
{
    return (poly_len - 2) * TRI_NPTS;
}

// https://en.wikipedia.org/wiki/Fan_triangulation
// keeps vertex orientation, returns triangles unchanged.
void fan_triangulate(vec3 *poly_pts, int poly_len, vec3 out_pts[fan_triangulate_out_len(poly_len)])
{
    int out_idx = 0;
    for (int i = 1; i < poly_len - 1; i++)
    {
        memcpy(out_pts[out_idx++], poly_pts[0], sizeof(vec3));

        // combining two memcpy calls
        memcpy(out_pts[out_idx], poly_pts[i], sizeof(vec3) * 2);
        out_idx += 2;
    }
    /* assert(out_idx == fan_triangulate_len(poly_len)); */
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
void surface_normal(vec3 dst, vec3 triangle[TRI_NPTS])
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

