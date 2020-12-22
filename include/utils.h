#ifndef UTILS_H
#define UTILS_H

#include "linmath_d.h"

double inc(double i, double lim);
int clamp(int val, int min, int max);
double normalize(double val, double upper, double lower);
double gamma_correct(double val, double g);
double distance(vec3 a, vec3 b);
void mesh_bounds(vec3 *verts, size_t n_verts, vec3 out_min, vec3 out_max);
void mesh_centroid(vec3 dst, vec3 *verts, size_t n_verts);
void surface_normal(vec3 dst, vec3 *triangle);
void normalize_coords(vec3 *out_arr, vec3 *verts, size_t n_verts);
void line_midpoint(vec2 dst, vec2 a, vec2 b);
void subdivide_line(vec2 *dst, size_t out_len, vec2 line[2], int n_subdivs);
int coord_to_idx(int x, int y, int rows, int cols);

#endif
