#include "fixed_point_32.h"
#include <stdint.h>
#include <math.h>

fixed_t fixed_from_dbl(double d)
{
    return (fixed_t)round((d * FIXED_SCALE_FACTOR));
}

double fixed_to_dbl(fixed_t f)
{
    return (double)f / FIXED_SCALE_FACTOR;
}

fixed_t fixed_mult(fixed_t f1, fixed_t f2)
{
    return ((fixed_extend_t)f1 * (fixed_extend_t)f2) / FIXED_SCALE_FACTOR;
}

fixed_t fixed_divide(fixed_t f1, fixed_t f2)
{
    return ((fixed_extend_t)f1 * FIXED_SCALE_FACTOR) / f2;
}
