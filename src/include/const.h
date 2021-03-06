#ifndef CONST_H
#define CONST_H

#include <gsl/gsl_sf_zeta.h>
#include <gsl/gsl_sf_psi.h>
#include "input.h"

////////////////////////////////////////////
// 					Constants			  //
////////////////////////////////////////////

namespace constants{
	const double pi = 3.14159265359 ;
	const double CA = 3.0, CF = 4.0/3, TF = 0.5, NC = 3.; 
	const double zeta3 = gsl_sf_zeta_int(3), zeta5 = gsl_sf_zeta_int(5) ;
	const double Euler = - gsl_sf_psi(1.0) ;
	const double pbarn = 3.894e8 ;
}

// SM paraneters
namespace smPara {
	extern const int Nf;
	extern const double MW, MW2, MZ, MZ2 ;
	extern const double SW, SW2, CW, CW2 ;
	extern const double GF, alpha, alpha2 ;
}

// Beam parameters
namespace beamPara{
	extern const double ECM, R, pTveto ;
}

// RG coefficients
namespace RGpara {
	extern const double B0, B1, B2; 
	extern const double G0, G1, G2; 
	extern const double gV0, gV1, gV2; 
	extern const double cLA, c0A, c2A, c4A; 
	extern const double cLf, c0f, c2f, c4f; 
}

////////////////////////////////////////////
#endif
