/* Try to avoid macros (especially global ones). If you can't do without
 * dump them here */

#if __STDC_VERSION__ < 199901L
# error Recompile with C99 support
#endif

#define Assert(...) Assert_Info(__func__, __FILE__, __LINE__, __VA_ARGS__)
#define Reallocate_P(...) Reallocate_P_Info(__func__, __FILE__, __LINE__, __VA_ARGS__)

#define Malloc(x) Malloc_info( __func__, __FILE__, __LINE__, x)
#define Realloc(x,y) Realloc_info(__func__, __FILE__, __LINE__, x, y)
#define Free(x) Free_info(__func__, __FILE__, __LINE__, x)
#define Profile(x) Profile_Info(__func__, __FILE__, __LINE__, x)

#define min(a,b) ((a)<(b)?(a):(b)) // this doesnt always work: c = max(a++, b)
#define max(a,b) ((a)>(b)?(a):(b))

#define len3(a) sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]) // these are slow !
#define len2(a) sqrt(a[0]*a[0] + a[1]*a[1])

#define len3_sq(a) (a[0]*a[0] + a[1]*a[1] + a[2]*a[2])
#define len2_sq(a) (a[0]*a[0] + a[1]*a[1])

#define p2(a) ((a)*(a))
#define p3(a) ((a)*(a)*(a))


