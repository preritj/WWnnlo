#ifndef BEAM_H
#define BEAM_H

#include "flavor.h"
#include "input.h"
#include "dist.h"

struct BeamPara{
	double mu; 		 //factorization scale
	double M; 		 // invariant mass of lepton pair
	double tau; 	 // ratio : (M/ECM)
	char flavor; 	 // flavor type : either 'u' or 'd' 
};

struct zPara{
	double mu; 		 //factorization scale
	double tau; 	 // ratio : (M/ECM) 
	double Y; 		 // rapidity 
	char flavor; 	 // flavor type : either 'u' or 'd'
	int flag;		 // relevant for z integration 
};

class Beam{
	private : 
		BeamPara para;
		double as, xsNNLO ;
	public :
		Beam() ; 
		void set_mu(double);
		flav lo(double), nlo(double) ;
};

double integrand_rap_lo(double, void*);
double integrand_rap_nlo(double, void*);
double integrand_z(double, void*);

#endif
