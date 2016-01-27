#include <iostream> 
#include <cstdlib> 
#include <math.h>
#include <const.h>
#include <input.h>
#include <PDF.h>
#include <Hard.h>
#include <misc.h>

using namespace std;
using namespace constants;

void Hard::set_mu(double mu_){
	mu=mu_; 
}

//----------------------------------------
// Debug functions 
flav Hard::H0(){
	return proc.Hlo(); 
}

flav Hard::H1reg(){
	return proc.Hnlo_reg(); 
}

//----------------------------------------
// LO Hard function 
flav Hard::lo(double M){
	flav result, error;
	hPara para(proc); 
	para.order= 0;
	para.M = M;

	//up-type 
	para.flavor = 'u';
	NIntegrate(0., 1., &para, integrand_cosT, result.u, error.u, 0);
	//up-type 
	para.flavor = 'd';
	NIntegrate(0., 1., &para, integrand_cosT, result.d, error.d, 0);

	return result; 
}

