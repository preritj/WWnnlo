#include <const.h>
#include <math.h>

//one-loop hard function
//---------------------------------------------------------
double H1(double M, double mu){
	double L = 2.*log(mu/M) ;
	//H1mu : log(mu/M) dependent pieces pieces
	double H1mu = -G0 * L*L/2 + 2.*g0*L ;  
	//H1reg: mu independent piece 
	return H1mu + H1reg(M) ;
}
//---------------------------------------------------------

//two-loop hard function

