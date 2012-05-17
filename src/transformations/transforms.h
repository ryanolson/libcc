///
// CC Tensor Index
///
typedef struct {
    int i[4];
    int n[4];
} ccti_t;


/// 
// External Subroutines
//
void translate_2d_to_4d(const ccti_t *t2d, ccti_t *t4d);
void translate_4d_to_2d(const ccti_t *t4d, ccti_t *t2d);
void increment_4d(ccti_t *t);
void increment_2d(ccti_t *t);

/// 
// Index Conversions
///
static inline swap(int i, int j)
{
        int tmp = i;
        i = j;
        j = tmp;
}

static inline void ccti_tran12(ccti_t *t) { swap(t->i[0], t->i[1]); }
static inline void ccti_tran13(ccti_t *t) { swap(t->i[0], t->i[2]); }
static inline void ccti_tran14(ccti_t *t) { swap(t->i[0], t->i[3]); }
static inline void ccti_tran23(ccti_t *t) { swap(t->i[1], t->i[2]); }
static inline void ccti_tran24(ccti_t *t) { swap(t->i[1], t->i[3]); }
static inline void ccti_tran34(ccti_t *t) { swap(t->i[2], t->i[3]); }

static inline void ccti_swap12(ccti_t *t) { ccti_tran12(t); swap(t->n[0], t->n[1]); }
static inline void ccti_swap13(ccti_t *t) { ccti_tran13(t); swap(t->n[0], t->n[2]); }
static inline void ccti_swap23(ccti_t *t) { ccti_tran23(t); swap(t->n[1], t->n[2]); }

static inline void ccti_insi12(ccti_t *t) { ccti_swap12(t); }
static inline void ccti_insi13(ccti_t *t) { ccti_swap13(t); }
static inline void ccti_insi14(ccti_t *t) { ccti_swap14(t); }

static inline void ccti_init4d(ccti_t *t, int i1, int i2, int i3, int i4, int n1, int n2, int n3, int n4)
{
        t->i[0] = i1;
        t->i[1] = i2;
        t->i[2] = i3;
        t->i[3] = i4;
        t->n[0] = n1;
        t->n[1] = n2;
        t->n[2] = n3;
        t->n[3] = n4;
}

static inline void ccti_init2d(ccti_t *t, int i1, int i2, int n1, int n2)
{
        ccti_init4d(t, i1, i2, 0, 0, n1, n2, 0, 0);
}

static inline void ccti_init1d(ccti_t *t, int i, int n)
{
        ccti_init2d(t, i, 0, n, 0);
}

static inline void ccti_get_4d_indexes(ccti_t *t, int *i1, int *i2, int *i3, int *i4)
{
        *i1 = t->i[0];
        *i2 = t->i[1];
        *i3 = t->i[2];
        *i4 = t->i[3];
}

static inline void ccti_get_4d_sizes(ccti_t *t, int *n1, int *n2, int *n3, int *n4)
{
        *n1 = t->n[0];
        *n2 = t->n[1];
        *n3 = t->n[2];
        *n4 = t->n[3];
}


static inline void ccti_get_2d_indexes(ccti_t *t, int *i1, int *i2)
{
        *i1 = t->i[0];
        *i2 = t->i[1];
}

static inline void ccti_get_2d_sizes(ccti_t *t, int *n1, int *n2)
{
        *n1 = t->n[0];
        *n2 = t->n[1];
}

static inline void ccti_get_2d_indexes(ccti_t *t, int *i)
{
        *i = t->i[0];
}

static inline void ccti_get_2d_sizes(ccti_t *t, int *n)
{
        *n = t->n[0];
}

