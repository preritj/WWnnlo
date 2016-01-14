#ifndef ANOM_DIM_H
#define ANOM_DIM_H

#include "const.h"

// constants that appear in beta function
const double B0 = (11*CA/3 - 4*TF*Nf/3) ;
const double B1 = (34*pow(CA,2)/3 - 20*CA*TF*Nf/3 - 4*CF*TF*Nf) ;
const double B2 = 2857*pow(CA,3)/54 
				+ (2*pow(CF,2) - 205*CF*CA/9 - 1415*pow(CA,2)/27)*
				TF*Nf + (44*CF/9 + 158*CA/27)*pow(TF*Nf,2) ;
// constants that appear in cusp anomalous dimension
const double G0 = (4*CF) ;
const double G1 = 4*CF*((67.0/9 - pow(pi,2)/3)*CA - 20.0/9*TF*Nf) ;
const double G2 = 4*CF*(pow(CA,2)*(245.0/6 -134*pow(pi,2)/27 +
					11*pow(pi,4)/45 + 22*zeta3/3) + CA*TF*Nf*
					(-418.0/27 + 40*pow(pi,2)/27 -56*zeta3/3) + 
					CF*TF*Nf*(-55.0/3+16*zeta3)-16*pow(TF*Nf,2)/27);

// constants that appear in quark anomalous dimension
const double gV0 = (-6*CF) ;
const double gV1 = pow(CF,2)*(-3.0 +4*pow(pi,2) -48*zeta3) + 
				+ CF*CA*(-961.0/27 -11*pow(pi,2)/3+52*zeta3)
				+ CF*TF*Nf*(260.0/27+4*pow(pi,2)/3);
const double gV2 = pow(CF,3)*(-29.0 - 6*pow(pi,2) - 16*pow(pi,4)/5 
					- 136*zeta3 + 32*pow(pi,2)/3*zeta3 + 480*zeta5 ) 
			+pow(CF,2)*CA*(-151.0/2 + 410*pow(pi,2)/9+494*pow(pi,4)/135
				 - 1688*zeta3/3 - 16*pow(pi,2)*zeta3/3 - 240*zeta5)
			+CF*pow(CA,2)*(-139345.0/1458 -7163*pow(pi,2)/243 -272*zeta5
			   - 83*pow(pi,4)/45 + 7052*zeta3/9 - 88*pow(pi,2)/9*zeta3)
			+pow(CF,2)*TF*Nf*(5906.0/27 -52*pow(pi,2)/9 -56*pow(pi,4)/27
				     + 1024*zeta3/9) 
			+CF*CA*TF*Nf*( -34636.0/729 + 5188*pow(pi,2)/243 
				      + 44*pow(pi,4)/45 - 3856*zeta3/27)
			+CF*pow(TF*Nf,2)*(19336.0/729-80*pow(pi,2)/27-64*zeta3/27) ;

// constants that appear in jet-veto functions
const double cLA = 131./72 - pi*pi/6 - 11./6*log(2.) ;
const double c0A = -805./216 + 11*pi*pi/72 + 35./18*log(2.)
					+11./6*pow(log(2.),2) + zeta3/2 ;
const double c2A = 1429./172800 + pi*pi/48 + 13./180*log(2.) ;
const double c4A = -0.0225794 ;
const double cLf = -23./36 + 2./3*log(2.) ;
const double c0f = 157./108 - pi*pi/18 - 8./9*log(2.) 
					- 2./3*pow(log(2.),2) ;
const double c2f = 3071./86400 - 7./360*log(2.) ;
const double c4f = -0.000442544 ;

#endif
