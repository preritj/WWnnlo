#include <iostream> 
#include <cstdlib> 
#include <math.h>
#include <string>
#include <vector>
#include <algorithm>
#include <const.h>
#include <PDF.h>
#include <Beam.h>
#include <input.h>

using namespace std;
using namespace beamPara;

Beam::Beam(){
	set_mu(pTveto); // default scale = pT veto
}

void Beam::set_mu(double mu_){
	para.mu=mu_;
}

flav Beam::lo(double M_, double mu_){
	set_mu(mu_); return lo(M_);
}
 
flav Beam::nlo(double M_, double mu_){
	set_mu(mu_); return nlo(M_);
} 
