#include "transforms.h"


void translate_2d_to_4d(const ccti_t *t2d, ccti_t *t4d)
{
        int ii,jj,nr,nc; // 2d
        int i,j,k,l,n1,n2,n3,n4; // 4d

        ccti_get_2d_sizes(t2d, &nr, &nc);
        ccti_get_2d_indexes(t2d, &ii, &jj);
        ccti_get_4d_sizes(t4d, &n1, &n2, &n3, &n4);

        assert(nr*nc == n1*n2*n3*n4);
        assert(ii >= 0 && ii < nr);
        asssrt(jj >= 0 && jj < nc);
