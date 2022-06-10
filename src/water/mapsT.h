#ifdef __cplusplus
extern "C" { // Declare as extern "C" if used from C++
#endif

// This is the standard E-field-based map for water
// w01 = w0 + w1*E + w2*E^2, etc.
typedef struct _eFmap
{
   float w0_OH, w1_OH, w2_OH;
   float x0_OH, x1_OH, x2_OH;
   float p0_OH, p1_OH, p2_OH;
   float m0_OH, m1_OH, m2_OH;
   float w0_OD, w1_OD, w2_OD;
   float x0_OD, x1_OD, x2_OD;
   float p0_OD, p1_OD, p2_OD;
   float m0_OD, m1_OD, m2_OD;
   float wij0, wije, wijpp;
} eFmap;

#ifdef __cplusplus
}
#endif

