#include <stdlib.h>
#include <stdio.h>
#include "lib/linmath_d.h"

int main(void)
{
    vec3 x = {1.,2.,3.};
    vec3 y = {3.,4.,5.};
    double z = vec3_mul_inner(x,y);

    printf("%f",z);
}
