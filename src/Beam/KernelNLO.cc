#include <iostream> 
#include <cstdlib> 
#include <math.h>
#include <string>
#include <vector>
#include <algorithm>
#include <const.h>
#include <functions.h>
#include <Beam.h>
#include <misc.h>

using namespace std;

Dist operator+(const Dist& d1, const Dist& d2){
	Dist d;
	d.None = d1.None + d2.None;
	d.Plus = d1.Plus + d2.Plus;
	d.Delta = d1.Delta + d2.Delta;
	return d;
}

Dist operator*(const Dist& d1, double factor){
	Dist d;
	d.None = d1.None*factor;
	d.Plus = d1.Plus*factor;
	d.Delta = d1.Delta*factor;
	return d;
}

Dist Iqfromq(double z, double Lp){
	Dist Pqfromq, Rqfromq ;
	Pqfromq.Plus = 4.*CF*(1+z*z)/(1-z) ;
	Rqfromq.None = CF*2.*(1-z) ;
	Rqfromq.Delta = -CF*pi*pi/6 ;
	return Pqfromq*(-Lp/2) + Rqfromq;
}


Dist Iqfromg(double z, double Lp){
	Dist Pqfromg, Rqfromg ;
	Pqfromg.None = 4.*TF*(z*z + (1-z)*(1-z)) ;
	Rqfromg.None = 4.*TF*z*(1-z) ;
	return Pqfromg*(-Lp/2) + Rqfromg;
}

double PlusRem(double xi){
	return 2*CF*(xi*(2+xi) + 4*log(1-xi)) ;
}

