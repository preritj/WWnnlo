#include <const.h>
#include <math.h>
#include "flav.h"

//one-loop hard function
//----------------------------------------------------------
flav H1(double M, double mu){
	double L = 2.*log(mu/M) ;
	//H1mu : log(mu/M) dependent pieces (flavor independent)
	double H1mu = -G0 * L*L/2 + 2.*g0*L ;  
	//H1reg: mu independent piece 
	flav _H1;
	flav _H1reg = H1reg(M);
	_H1.u = H1mu + _H1reg.u ;
	_H1.d = H1mu + _H1reg.d ;
	return _H1;
}

flav H1reg(double M){
	double c0 = pi*pi*

}
//----------------------------------------------------------

//two-loop hard function

