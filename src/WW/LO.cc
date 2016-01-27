#include <const.h>
#include <WW.h>
#include <math.h>

using namespace smPara;

double pp2WW::F0(){
	return
	(16*s)/MW2 + 16*(0.25 + pow(MW2,2)/pow(t,2))*(-1 + 
	(t*u)/pow(MW2,2));
}

double pp2WW::J0(){
	return
	16*(-2 + s/MW2 + (2*MW2)/t) + 16*(-MW2/s/2. + 1./4. - 
	pow(MW2,2)/t/s)*(-1 + (t*u)/pow(MW2,2));
}

double pp2WW::K0(){
	return
	8*(-4 + s/MW2) + 8*(3*pow(MW2/s,2) - MW2/s + 1./4.)*(-1 
	+ (t*u)/pow(MW2,2));
}

flav pp2WW::CQst(){
	return 
	-4*SW2*(Q + (T3-Q*SW2)*(s/(s-MZ2)/SW2))*PosUp;
}

flav pp2WW::CQss(){
	return 
	16*SW2*SW2*( pow(Q + (T3*0.5-Q*SW2)*(s/(s-MZ2)/SW2)  ,2) 
		+ pow( T3*(s*0.5/(s-MZ2)/SW2)  ,2.)  ) ;
}
