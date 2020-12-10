#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <stdint.h>
#include <assert.h>
#include <math.h>

// https://github.com/datenwolf/linmath.h
//   i made all the vec# types doubles, added string.h
#include "lib/linmath_d.h"

#include "obj_parser.h"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "lib/stb_image_write.h"

double normalize(double val, double upper, double lower);

int orient2d(vec2 a, vec2 b, vec2 c);
void triangle(vec3 a, vec3 b, vec3 c, int rows, int cols,
        uint8_t buffer[rows][cols], double z_buffer[rows][cols],
        double subpixel, uint8_t fill);

void surface_normal(vec3 dst, vec3 a, vec3 b, vec3 c);
void camera_transform(mat4x4 dst, vec3 camera_pos, vec3 target, vec3 up);

void dump_array(const char *tag, int rows, int cols, uint8_t array[rows][cols]);


int main(void)
{
    size_t size;
    double (*verts)[3];
    double (*texcoords)[2];
    double (*normals)[3];

    parse_obj("models/cube_quad.obj", &size, &verts, &texcoords, &normals);

    // parse_obj should return all faces as tris, but let's make sure
    assert(size % 3 == 0);

    int rows = 2048;
    int cols = 2048;
    uint8_t (*buffer)[cols] = calloc(rows * cols, sizeof(uint8_t));
    double (*z_buffer)[cols] = calloc(rows * cols, sizeof(double));

    mat4x4 m;
    vec3 up = {0, 1, 0};
    vec3 cpos = {0,0,0};
    camera_transform(m, cpos, verts[0], up);

    // loop through each face
    for (int i = 0; i < size; i+=3)
    {

        // world coords
        double *v1 = verts[i]; // x y z
        double *v2 = verts[i+1];
        double *v3 = verts[i+2];

        vec4 projp;
        mat4x4_mul_vec4(projp, m, v1);

        double up = 1;
        double low = -1;
        vec2 scrp = {
            normalize(projp[0], up, low),
            normalize(projp[1], up, low),
        };

        // word coords but homogeneous
        /* vec4 v1_h = {v1[0], v1[1], v1[2], 1}; */
        /* vec4 v2_h = {v2[0], v2[1], v2[2], 1}; */
        /* vec4 v3_h = {v3[0], v3[1], v3[2], 1}; */

        // homogeneous screen coords, packing z into 3rd space
        /* vec4 sp1 = {0, 0, 0, 1}; */
        /* vec4 sp2 = {0, 0, 0, 1}; */
        /* vec4 sp3 = {0, 0, 0, 1}; */
        /* sp1[2] = v1[2]; */
        /* sp2[2] = v2[2]; */
        /* sp3[2] = v3[2]; */
        /* orthographic_projection(v1_h, sp1); */
        /* orthographic_projection(v2_h, sp2); */
        /* orthographic_projection(v3_h, sp3); */


        // computing screen coords
        /* double p1_x = ((v1[0] + 1.0) * cols / 2.0); */
        /* double p1_y = ((v1[1] + 1.0) * rows / 2.0); */

        /* double p2_x = ((v2[0] + 1.0) * cols / 2.0); */
        /* double p2_y = ((v2[1] + 1.0) * rows / 2.0); */

        /* double p3_x = ((v3[0] + 1.0) * cols / 2.0); */
        /* double p3_y = ((v3[1] + 1.0) * rows / 2.0); */

        // storing z value of original point in screen space point
        /* vec3 screenpts[3]; */
        /* screenpts[0][0] = p1_x; */
        /* screenpts[0][1] = p1_y; */
        /* screenpts[0][2] = v1[2]; */

        /* screenpts[1][0] = p2_x; */
        /* screenpts[1][1] = p2_y; */
        /* screenpts[1][2] = v2[2]; */

        /* screenpts[2][0] = p3_x; */
        /* screenpts[2][1] = p3_y; */
        /* screenpts[2][2] = v3[2]; */

        /* vec3 normal; */
        /* surface_normal(normal, v1, v2, v3); */
        /* vec3_norm(normal, normal); */

        /* vec3 light_direction = {0, 0, -1}; // negative camera direction */
        /* float intensity = vec3_mul_inner(light_direction, normal); */

        // backface culling
        /* if (intensity > 0) */
        /* { */
        /*     triangle(screenpts[0], screenpts[1], screenpts[2], */
        /*             rows, cols, buffer, z_buffer, 0.5, intensity * 255); */
        /* } */
    }

    stbi_write_jpg("./fuck.jpg",rows, cols, 1, buffer, 100);

    free(verts);
    free(texcoords);
    free(normals);
    free(buffer);
    free(z_buffer);

}

double normalize(double val, double upper, double lower)
{
    return (val - lower) / (upper - lower);
}

//https://stackoverflow.com/questions/1560492/how-to-tell-whether-a-point-is-to-the-right-or-left-side-of-a-line
int orient2d(vec2 a, vec2 b, vec2 c)
{
    return (b[0]-a[0])*(c[1]-a[1]) - (b[1]-a[1])*(c[0]-a[0]);
}


// find the bbox of a triangle, loop over each pixel in bbox,
//    see if it's in the triangle by checking for negative bary coords
void triangle(vec3 a, vec3 b, vec3 c, int rows, int cols,
        uint8_t buffer[rows][cols], double z_buffer[rows][cols],
        double subpixel, uint8_t fill)
{
    // getting bounding box of the triangle
    double max_x = fmax(a[0], fmax(b[0], c[0]));
    double max_y = fmax(a[1], fmax(b[1], c[1]));
    double min_x = fmin(a[0], fmin(b[0], c[0]));
    double min_y = fmin(a[1], fmin(b[1], c[1]));

    vec3 sp; // screen space point
    // looping through every point in the bounding box
    for (sp[0] = min_x; sp[0] <= max_x; sp[0] += subpixel)
    {
        for (sp[1] = min_y; sp[1] <= max_y; sp[1] += subpixel)
        {
            // if bbox point is out of screen bounds
            if (sp[0] >= cols - 1 || sp[1] >= rows - 1 || sp[0] < 0 || sp[1] < 0)
            {
                continue;
            }

            // barycentric coords
            int w0 = orient2d(b, c, sp);
            int w1 = orient2d(c, a, sp);
            int w2 = orient2d(a, b, sp);

            // if point is inside or on all edges
            if (w0 >= 0 && w1 >= 0 && w2 >= 0)
            {
                sp[2] = 0; // init z value for the z buffer

                // smearing actual z values over bayercentric coords
                sp[2] += a[2] * w0;
                sp[2] += b[2] * w1;
                sp[2] += c[2] * w2;

                // flipping vertically
                int row = (rows - 1) - ceil(sp[1]);
                int col = ceil(sp[0]);

                if (z_buffer[row][col] < sp[2])
                {
                    z_buffer[row][col] = sp[2];
                    buffer[row][col] = fill;
                }
            }
        }
    }
}

// computes surface normal of triangle abc
void surface_normal(vec3 dst, vec3 a, vec3 b, vec3 c)
{
    vec3 v;
    vec3 w;
    vec3_sub(v, c, a);
    vec3_sub(w, b, a);
    vec3_mul_cross(dst, v, w);
}

// https://www.3dgep.com/understanding-the-view-matrix/#Transformations
void camera_transform(mat4x4 dst, vec3 camera_pos, vec3 target, vec3 up)
{
    vec3 z_axis; // camera forward vector
    vec3_sub(z_axis, camera_pos, target);
    vec3_norm(z_axis, z_axis);

    vec3 x_axis; // camera right vector
    vec3_mul_cross(x_axis, up, z_axis);
    vec3_norm(x_axis, x_axis);

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

void dump_array(const char *tag, int rows, int cols, uint8_t array[rows][cols])
{
    printf("%s (%dx%d):\n", tag, rows, cols);
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
        {
            printf("%3d", array[i][j]);
        }
        putchar('\n');
    }
}



