#include <WW.h>

using namespace constants;
using namespace smPara;

pp2WW::pp2WW(){
	M1=MW; M2=MW;
	Q.u = 2./3 ; Q.d = -1./3;
	T3.u = 1./2 ; T3.d = -1./2 ;
	PosUp.u = 1.; PosUp.d = -1. ;
	coup = pi*pi*alpha2 ;
	CQtt.u = 1. ; CQtt.d = 1. ;
}

void pp2WW::calc_lo(){
	H0 = 1./4/NC*coup*(CQtt*F0() + CQss()*K0() + CQst()*J0());
}

void pp2WW::calc_nlo(){
	H1 = CF/8./NC*coup*(CQtt*F1() + CQss()*K1() + CQst()*J1());
}
