#ifndef HARD_H
#define HARD_H

#include "flavor.h"

class Hard{
	private : 
		double MW, MZ, MW2, MZ2;
		double SW, CW, SW2, CW2;
		double GF, alpha, alpha2;
		double as, xsNNLO ;
	public : 
		Hard(InputPara);
		void set_mu(double);
		flav lo(double), nlo(double) ;
};
#endif
