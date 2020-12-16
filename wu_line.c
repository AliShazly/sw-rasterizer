#include <math.h>
#include <stdint.h>
#include <stdbool.h>
#include <assert.h>

// Implimentation of Xiaolin Wu's line algorithm
// https://en.wikipedia.org/wiki/Xiaolin_Wu%27s_line_algorithm

void swap(void *a, void *b, int size)
{
    char tmp; // single byte buffer
    for (int i = 0; i < size; i++)
    {
        tmp = *((char*)a + i);
        *((char*)a + i) = *((char*)b + i);
        *((char*)b + i) = tmp;
    }
}

void plot(int x, int y, double c, int rows, int cols, uint8_t (*out_buf)[3], uint8_t color[3])
{
    // out of bounds
    if ((x > cols - 1 || x < 0) || (y > rows - 1 || y < 0))
    {
        return;
    }

    // flipping vertically
    int row = (rows - 1) - y;
    int col = x;

    // 2D index (row,col) to a 1D index
    int idx = row * cols + col;
    for (int i = 0; i < 3; i++)
    {
        out_buf[idx][i] = c * color[i];
    }
}

double ipart(double x)
{
    return floor(x);
}

double round_(double x)
{
    return ipart(x + 0.5);
}

double fpart(double x)
{
    return x - floor(x);
}

double rfpart(double x)
{
    return 1 - fpart(x);
}

void draw_line(double x0, double y0, double x1, double y1,
        int rows, int cols, uint8_t (*out_buf)[3], uint8_t color[3])
{
    bool steep = fabs(y1 - y0) > fabs(x1 - x0);
    if (steep)
    {
        swap(&x0, &y0, sizeof(double));
        swap(&x1, &y1, sizeof(double));
    }
    if (x0 > x1)
    {
        swap(&x0, &x1, sizeof(double));
        swap(&y0, &y1, sizeof(double));
    }

    double dx = x1 - x0;
    double dy = y1 - y0;
    double gradient = dy / dx;
    if (dx == 0.)
    {
        gradient = 1.;
    }

    double xend = round_(x0);
    double yend = y0 + gradient * (xend - x0);
    double xgap = rfpart(x0 + 0.5);
    double xpxl1 = xend;
    double ypxl1 = ipart(yend);
    if (steep)
    {
        plot(ypxl1,   xpxl1, rfpart(yend) * xgap, rows, cols, out_buf, color);
        plot(ypxl1+1, xpxl1,  fpart(yend) * xgap, rows, cols, out_buf, color);
    }
    else
    {
        plot(xpxl1, ypxl1  , rfpart(yend) * xgap, rows, cols, out_buf, color);
        plot(xpxl1, ypxl1+1,  fpart(yend) * xgap, rows, cols, out_buf, color);
    }
    double intery = yend + gradient;

    xend = round_(x1);
    yend = y1 + gradient * (xend - x1);
    xgap = fpart(x1 + 0.5);
    double xpxl2 = xend;
    double ypxl2 = ipart(yend);
    if (steep)
    {
        plot(ypxl2  , xpxl2, rfpart(yend) * xgap, rows, cols, out_buf, color);
        plot(ypxl2+1, xpxl2,  fpart(yend) * xgap, rows, cols, out_buf, color);
    }
    else
    {
        plot(xpxl2, ypxl2,  rfpart(yend) * xgap, rows, cols, out_buf, color);
        plot(xpxl2, ypxl2+1, fpart(yend) * xgap, rows, cols, out_buf, color);
    }

    if (steep)
    {
        for (int x = xpxl1 + 1; x < xpxl2; x++)
        {
            plot(ipart(intery)  , x, rfpart(intery), rows, cols, out_buf, color);
            plot(ipart(intery)+1, x,  fpart(intery), rows, cols, out_buf, color);
            intery = intery + gradient;
        }
    }
    else
    {
        for (int x = xpxl1 + 1; x < xpxl2; x++)
        {
            plot(x, ipart(intery),  rfpart(intery), rows, cols, out_buf, color);
            plot(x, ipart(intery)+1, fpart(intery), rows, cols, out_buf, color);
            intery = intery + gradient;
        }
    }
}

