#include <stdio.h>
#include "cuda_runtime_api.h"

#define INDX(a,b,c,ld) ( ( (c) * (ld) * (ld) ) \
                       + ( (b) * (ld) ) \
                       + ( (a) ) )

#define SHARED_REDUCTION_SIZE 128
#define NOTEX


texture<int2,1>  tex_x_double_v3;

static __inline__ __device__ double fetch_x_v3(const int& i)
{
  register int2  v = tex1Dfetch(tex_x_double_v3, i);
  return __hiloint2double(v.y, v.x);
}
