#ifndef HARD_H
#define HARD_H

#include "flavor.h"

void anom_dim(int) ;

class Hard{
	protected : 
		double as, mu ;
	public : 
		void set_mu(double);
		flav lo(), nlo() ;
		flav lo(double), nlo(double) ;
};
#endif
