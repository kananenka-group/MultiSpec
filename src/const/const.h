#ifndef CONST_H
#define CONST_H

#include <cmath>

namespace constants
{
// common constants
const double HBAR = 5.308837367; // in cm-1 * ps
const double AU_TO_WN = 2.19474e5;
const double A0 = 1.889726125;

// masses of common atoms
const float mO = 16.000;
const float mH = 1.008;

// WATER MODULE------------
//7.831 A * 1.88 conversion to au, 7.831 A cut-off is chosen for historical reasons
// see: JCTC 9, 3109 (2013), JCP 128, 224511 (2008), JCP 132, 204505 (2010)
const float O_to_H_dist_water_cutoff = 14.798445;

// Ratio between the longitudinal and transverse bond polarizability derivatives
// see J. G. Scherer and R. G. Snyder, J. Chem. Phys. 67, 4794 (1977).
const float bond_plz_ratio = 5.6;

// Equilibrium constant for: H2O + D2O = 2HOD
// K = 3.828 is taken as an average from:
// Table 1 of:
// Yongdan Kim, Namkung Hankyu, and Hoeil Chung,
// "Determination of the Equilibrium Constant of the Isotropic Disproportionation
//  Between Water and Heavy Water Using Near-Infrared Spectroscopy,"
//  Appl. Spectrosc. 63, 256-259 (2009)
//#define  kEqIsoW  3.9
// Another possible value is 3.9 from:
// ``Orthogonalyzed H2O and D2O species obtained from infrared spectra of liquid water at several temperatures''
// by Jean-Joseph Max, Pascal Larouche and Camille Chapados
// Journal of Molecular Structure 1149 (2017) 457-472
const float kEqIsoW = 3.9;

// Switching function cut-off for SFG calculations (in A)
// See J. Chem. Phys. 135, 044701 (2011)
const float SWITCHF_CUT = 4.0;

// AMIDE I MODULE----------------
const float CN_dist_cutoff = 20.0*A0;      //20A
const float amideItdAngle = 10.0*(M_PI/180.0);
const float amideItdMag = 2.73;  //D/A/amu^1/2
const float amideItd1 = 1.256668374; // 0.665 A -> a0
const float amideItd2 = 0.487549535; // 0.258 A -> a0
const float amideIcoupConst = 383.313*A0*A0*A0/(amideItdMag*amideItdMag);
}


#endif
