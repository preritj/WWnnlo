#include <iostream> 
#include <cstdlib> 
#include <math.h>
#include <string>
#include <vector>
#include <algorithm>
#include <const.h>
#include <PDF.h>
#include <const.h>
#include <Beam.h>
#include <misc.h>

using namespace std;
using namespace smPara;
using namespace beamPara;

//-----------------------------------------------------------------
// LO product of Beam functions
flav Beam::lo(double M_){
	para.M = M_;
	para.tau = M_/ECM ;
	flav result, error;

	// up-type 
	para.flavor = 'u';
	NIntegrate(0., 1., &para, integrand_rap_lo, result.u, error.u, 0);
	//down-type 
	para.flavor = 'd';
	NIntegrate(0., 1., &para, integrand_rap_lo, result.d, error.d, 0);
	return result;
}

//-----------------------------------------------------------------
//  Integrand for the rapidity integral (limits from 0 to 1)
double integrand_rap_lo(double x, void* para_ptr)
{
	double Jac = 1.0 ;  // Take care of Jacobian in the integrals
	BeamPara para = *(BeamPara*) para_ptr ;
	double tau = para.tau; double mu = para.mu ;

	vector<int> flav;
	if (para.flavor == 'u') {flav.push_back(2); flav.push_back(4);} 
    else if (para.flavor == 'd') {
		flav.push_back(1); flav.push_back(3);
		if (Nf == 5) flav.push_back(5); }
	else {cout << "Wrong flavor type. Aborting..." << endl; abort();}
	
	double Y_min = log(tau) ;
	double Y_max =  -log(tau) ;
	double Y = Y_min + x*(Y_max - Y_min) ;
	Jac = Jac * (Y_max - Y_min) ;
	
	// Call PDFs (note that function returns x*f(x) )
	double f1[13], f2[13] ;
	evolvePDF(tau*exp(Y), mu, f1) ; 
	evolvePDF(tau*exp(-Y), mu, f2) ;

	double integrand = 0. ;
	for (int i=0; i < flav.size(); i++){
		integrand += f1[6+ flav[i] ] * f2[6 -flav[i] ] /tau/tau ; 
	}
	return Jac*integrand ;
}

