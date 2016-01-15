#ifndef CONST_H
#define CONST_H

#include <gsl/gsl_sf_zeta.h>
#include <gsl/gsl_sf_psi.h>
#include "input.h"

////////////////////////////////////////////
// 					Constants			  //
////////////////////////////////////////////
const double pi = 3.14159265359 ;
const double CA = 3.0, CF = 4.0/3, TF = 0.5, NC = 3.; 
const double zeta3 = gsl_sf_zeta_int(3), zeta5 = gsl_sf_zeta_int(5) ;
const double Euler = - gsl_sf_psi(1.0) ;
const double pbarn = 3.894e8 ;

// derived constants from input parameters
const double alpha = 1.0/130 ;
//double MW2 = MW*MW; 
//const double MW2 = MW*MW ;
//const double MZ2 = MZ*MZ ;
//const double SW2 = 1.-MW2/MZ2, SW = sqrt(SW2) ;
//const double CW2 = 1.-SW2, CW = sqrt(CW2) ;

////////////////////////////////////////////
#endif
