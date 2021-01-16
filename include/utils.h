#ifndef UTILS_H
#define UTILS_H

#include "linmath_d.h"
#include <stdbool.h>

/* 3 points in a triangle */
#define TRI_NPTS 3

void print_vecn(double* vecn, int n);
double inc(double i, double lim);
int ptr_dist(void *p1, void* p2);
int clamp(int val, int min, int max);
int imin(int a, int b);
int imax(int a, int b);
bool dbl_almost_equal(double a, double b, double eps);
double normalize(double val, double upper, double lower);
double gamma_correct(double val, double g);
double distance(vec3 a, vec3 b);
double lerp(double a, double b, double f);
double inv_lerp(double a, double b, double x);
bool is_inbetween(double bound1, double bound2, double val);
bool point_in_bbox(vec2 p, double max_x, double min_x, double max_y, double min_y);
void plane_line_intersection(vec3 dst, vec3 a, vec3 b, vec3 plane_pt, vec3 plane_normal);
void mesh_bounds(vec3 *verts, size_t n_verts, vec3 out_min, vec3 out_max);
void mesh_centroid(vec3 dst, int n_verts, int vert_size, double (*verts)[vert_size]);
void triangle_normal(vec3 dst, int tri_pt_size, double triangle[TRI_NPTS][tri_pt_size]);
void normalize_coords(vec3 *out_arr, vec3 *verts, size_t n_verts);
void line_midpoint(vec2 dst, vec2 a, vec2 b);
void subdivide_line(vec2 *dst, size_t out_len, vec2 line[2], int n_subdivs);

#endif
