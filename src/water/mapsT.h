#ifdef __cplusplus
extern "C" { // Declare as extern "C" if used from C++
#endif

// This is the standard E-field-based map for water
// w01 = w0 + w1*E + w2*E^2, etc.
typedef struct _eFmap
{
   float w0, w1, w2;
   float x0, x1, x2;
   float p0, p1, p2;
   float m0, m1, m2;
   float wij0, wije, wijpp;
} eFmap;

#ifdef __cplusplus
}
#endif

