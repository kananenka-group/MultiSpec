#ifdef __cplusplus
extern "C" { // Declare as extern "C" if used from C++
#endif

//----------------------------------------------------
// This is the standard E-field-based map for water
// hydroxyl stretch
// ---------------------------------------------------
// w01 = w0 + w1*E + w2*E^2
// x01 = x0 + x1*w01 + x2*w01*w01
// p01 = p1 + p1*w01 + p2*w01*w01
// m01 = m0 + m1*E + m2*E^2
// Intramolecular coupling:
// wc  = [ wij0 +  wije*(Ei + Ej) ]xi*xj + wijpp*pi*pj
typedef struct _eFmapS
{
   float w10_OH_0, w10_OH_1, w10_OH_2;
   float x10_OH_0, x10_OH_1, x10_OH_2;
   float p10_OH_0, p10_OH_1, p10_OH_2;
   float m10_OH_0, m10_OH_1, m10_OH_2;
   float w10_OD_0, w10_OD_1, w10_OD_2;
   float x10_OD_0, x10_OD_1, x10_OD_2;
   float p10_OD_0, p10_OD_1, p10_OD_2;
   float m10_OD_0, m10_OD_1, m10_OD_2;
   float w21_OH_0, w21_OH_1, w21_OH_2;
   float x21_OH_0, x21_OH_1, x21_OH_2;
   float p21_OH_0, p21_OH_1, p21_OH_2;
   float m21_OH_0, m21_OH_1, m21_OH_2;
   float w21_OD_0, w21_OD_1, w21_OD_2;
   float x21_OD_0, x21_OD_1, x21_OD_2;
   float p21_OD_0, p21_OD_1, p21_OD_2;
   float m21_OD_0, m21_OD_1, m21_OD_2;
   float wij0, wije, wijpp;
} eFmapS;

#ifdef __cplusplus
}
#endif

