#ifndef __CONST_H__    //  __CONST_H__
#define __CONST_H__

// basic physical constants

#define pi  3.141592653589793238462643383279502884197169399375105
#define ev2eh 2.4179894887970428E14  // what's this in fundamental constants?

#if 1 // CGS
#define c   2.99792458E10            // cm/s
#define h   6.6260687652E-27         // erg s
#define k   1.3806503E-16            // erg/K
#define me  9.1093818872E-28         // g
#define mp  1.67262158E-24           // g
#define mn  1.67492716E-24           // g
#define mu  1.66053892173E-24        // g

#define e   4.8032042510E-10         // statC
//#define e   4.803204198984205264E-10 // statC? slight inconsistency here...
#define Ryd 1.097373156853955E5      // cm-1
#define e0  7.957747154594775607E-02 // 1/(4 pi)
#define R   8.31447215E7             // erg mol-1 K-1
#define ev  1.602189E-12             // erg
#define saha_const 2.0706615382E-16 // In cgs units

#define delta_T 1.0 // step for the computation of numerical derivatives
#define delta_vt 20.0 // step for the computation of numerical derivative
#define delta_Nt_frac 1E-3 
#define delta_vr 20.0 // same as for delta_vt
#define delta_B 0.01
#define delta_angle 0.0001 // In radians, this is a bit more than 0.1 degree. Hopefully it works. 

#define n_air = 1.000277 // in the end we took the VALC approach, so n is still lambda dependent. 

#define turn_on_damping 1.0 //

#else

#define c   2.99792458E8             // m/s
#define h   6.6260687652E-34         // J s
#define k   1.3806503E-23            // J/K
#define me  9.1093818872E-31         // kg
#define mp  1.67262158E-27           // kg
#define mn  1.67492716E-27           // kg
#define mu  1.66053892173E-27        // kg
#define e   1.60217646263E-19        // C
#define Ryd 1.097373156853955E7      // m-1
#define e0  8.854187817E-12          // F m-1
#define R   8.31447215               // J mol-1 K-1
#define ev  e                        // J

#endif

#endif                 //  __CONST_H__
