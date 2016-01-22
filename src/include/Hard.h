#ifndef HARD_H
#define HARD_H

#include "flavor.h"

void anom_dim(int) ;

class Hard{
	protected : 
		double MW, MZ, MW2, MZ2 ;
		double SW, SW2, CW, CW2 ;
		double GF, alpha, alpha2 ;
		double as, mu, Nf ;
	public : 
		Hard(InputPara);
		void set_mu(double);
		flav lo(), nlo() ;
		flav lo(double), nlo(double) ;
};
#endif
