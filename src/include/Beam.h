#ifndef BEAM_H
#define BEAM_H

#include "flavor.h"
#include "input.h"
#include "dist.h"

struct BeamPara{
	double ECM;		 // CoM energy
	double pTveto;   // pT-veto
	double mu; 		 //factorization scale
	double R; 		 // jet radius
	double M; 		 // invariant mass of lepton pair
	double tau; 	 // ratio : (M/ECM) 
	int Nf; 		 // number of flavors 
	char flavor; 	 // flavor type : either 'u' or 'd' 
};

struct zPara{
	double pTveto;   // pT-veto
	double mu; 		 //factorization scale
	double tau; 	 // ratio : (M/ECM) 
	double Y; 		 // rapidity 
	int Nf; 		 // number of flavors 
	char flavor; 	 // flavor type : either 'u' or 'd'
	int flag;		 // relevant for z integration 
};

class Beam{
	private : 
		BeamPara para;
		double as, xsNNLO ;
	public : 
		Beam(InputPara);
		void set_mu(double);
		flav lo(double), nlo(double) ;
};

double integrand_rap_lo(double, void*);
double integrand_rap_nlo(double, void*);
double integrand_z(double, void*);

#endif
