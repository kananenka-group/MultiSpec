#ifdef __cplusplus
extern "C" { // Declare as extern "C" if used from C++
#endif

// This is the standard E-field-based map for water
// bending vibration
typedef struct _eFmapB
{
   float w0_HOH_10, w1_HOH_10; 
   float w0_HOH_21, w1_HOH_21;
   float t0_HOH_10, t1_HOH_10;
   float t0_HOH_21, t1_HOH_21;
   float mup_HOH;
   float a_HOH_xx, a_HOH_yy, a_HOH_zz;

   float w0_DOD_10, w1_DOD_10; 
   float w0_DOD_21, w1_DOD_21;
   float t0_DOD_10, t1_DOD_10;
   float t0_DOD_21, t1_DOD_21;
   float mup_DOD;
   float a_DOD_xx, a_DOD_yy, a_DOD_zz;

   float w0_HOD_10, w1_HOD_10; 
   float w0_HOD_21, w1_HOD_21;
   float t0_HOD_10, t1_HOD_10;
   float t0_HOD_21, t1_HOD_21;
   float mup_HOD;
   float a_HOD_xx, a_HOD_yy, a_HOD_zz;
} eFmapB;

#ifdef __cplusplus
}
#endif

