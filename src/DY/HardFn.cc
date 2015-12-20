#include <const.h>
#include <math.h>

//one-loop hard function
//---------------------------------------------------------
double H1(double M, double mu){
	return CF/pi/2*(-8+ 7*pi*pi/6 + log(M/mu)*(6-4*log(M/mu))) ;
}
//---------------------------------------------------------

//two-loop hard function

