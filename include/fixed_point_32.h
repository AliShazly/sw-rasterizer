#ifndef FIXED_POINT_32_H
#define FIXED_POINT_32_H

#include <stdint.h>

// 28.4 fixed point
#define FIXED_FRACTION_BITS 4
#define FIXED_SCALE_FACTOR (1 << FIXED_FRACTION_BITS)

typedef int32_t fixed_t;
typedef int64_t fixed_extend_t;

fixed_t fixed_from_dbl(double d);
double fixed_to_dbl(fixed_t f);
fixed_t fixed_mult(fixed_t f1, fixed_t f2);
fixed_t fixed_divide(fixed_t f1, fixed_t f2);

#endif
