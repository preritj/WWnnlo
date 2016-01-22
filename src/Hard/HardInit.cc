#include <iostream> 
#include <cstdlib> 
#include <math.h>
#include <const.h>
#include <input.h>
#include <PDF.h>
#include <Hard.h>

using namespace std;

Hard::Hard(InputPara input){
	MW = input.MW; 
	MZ = input.MZ; 
	GF = input.GF;
	MW2=MW*MW; MZ2=MZ*MZ;
	SW2=1.-MW2/MZ2,		SW=sqrt(SW2);
	CW2=1.-SW2;			CW=sqrt(CW2);
	alpha = sqrt(2.)*MW2/pi*GF;
	alpha2 = alpha*alpha ;
	Nf = input.Nf ;
	anom_dim(Nf) ;
}

void Hard::set_mu(double mu_){
	mu=mu_; as=alpha_s(mu);
}

flav Hard::lo(double mu_){
	set_mu(mu_); return lo(); 
}

flav Hard::lo(){
	flav H; H.u = 1.; H.d = 1.; return H; 
}
