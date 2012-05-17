///
// CC Tensor Index
///
typedef struct {
    int i[4];
    int n[4];
} ccti_t;


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


