#include <const.h>
#include <WW.h>

using namespace smPara;

double WW::F0(){
	return
	(16*s)/MW2 + 16*(0.25 + (MW2*MW2)/(t*t))*(-1 + (t*u)/(MW2*MW2)) ;
}

double WW::J0(){
	return
	16*s*(-2 + s/MW2 + (2*MW2)/t) + 16*(-MW2/2. + s/4. 
	- (MW2*MW2)/t)*(-1 + (t*u)/(MW2*MW2)) ;
}

double WW::K0(){
	return
	8*(s*s)*(-4 + s/MW2) + 8*(3*(MW2*MW2) - MW2*s 
	+ (s*s)/4.)*(-1 + (t*u)/(MW2*MW2)) ;
}
