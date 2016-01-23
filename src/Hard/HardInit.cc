#include <iostream> 
#include <cstdlib> 
#include <math.h>
#include <const.h>
#include <input.h>
#include <PDF.h>
#include <Hard.h>

using namespace std;

void Hard::set_mu(double mu_){
	mu=mu_; as=alpha_s(mu);
}

flav Hard::lo(double mu_){
	set_mu(mu_); return lo(); 
}

flav Hard::lo(){
	flav H; H.u = 1.; H.d = 1.; return H; 
}
