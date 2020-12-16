#ifndef WU_LINE_H
#define WU_LINE_H

#include <stdint.h>

void swap(void *a, void *b, int size);
void plot(int x, int y, double c, int rows, int cols, uint8_t (*out_buf)[3], uint8_t color[3]);
double ipart(double x);
double round_(double x);
double fpart(double x);
double rfpart(double x);
void draw_line(double x0, double y0, double x1, double y1,
        int rows, int cols, uint8_t (*out_buf)[3], uint8_t color[3]);

#endif
