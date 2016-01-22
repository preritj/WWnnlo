#include <iostream> 
#include <cstdlib> 
#include <math.h>
#include <string>
#include <vector>
#include <algorithm>
#include <const.h>
#include <PDF.h>
#include <Beam.h>
#include <misc.h>

using namespace std;

//-----------------------------------------------------------------
// NLO product of Beam functions
flav Beam::nlo(double M_in){
	para.M = M_in;
	para.tau = M_in/para.ECM ;
	flav result, error;

	// up-type 
	para.flavor = 'u';
	NIntegrate(0., 1., &para, integrand_rap_nlo, result.u, error.u, 0);
	//down-type 
	para.flavor = 'd';
	NIntegrate(0., 1., &para, integrand_rap_nlo, result.d, error.d, 0);
	return result;
}

//-------------------------------------------------------------------
//  Integrand for the rapidity integral (limits from 0 to 1)
double integrand_rap_nlo(double x, void* para_ptr)
{
	double Jac = 1.0 ;  // Take care of Jacobian in the integrals
	BeamPara para = *(BeamPara*) para_ptr ;
	double tau = para.tau; double mu = para.mu ; double pTveto=para.pTveto;
	double Lp = 2*log(mu/pTveto) ;

	vector<int> flav; flav.clear();
	if (para.flavor == 'u') {flav.push_back(2); flav.push_back(4);} 
    else if (para.flavor == 'd') {
		flav.push_back(1); flav.push_back(3);
		if (para.Nf == 5) flav.push_back(5); }
	else {cout << "Wrong flavor type. Aborting..." << endl; abort();}
	
	double Y_min = log(tau) ;
	double Y_max =  -log(tau) ;
	double Y = Y_min + x*(Y_max - Y_min) ;
	Jac = Jac * (Y_max - Y_min) ;

	// Call PDFs (note that function returns x*f(x) )
	double f1[13], f2[13] ;
	double x1 = tau*exp(Y) ;
	double x2 = tau*exp(-Y) ;
	evolvePDF(x1, mu, f1) ; 
	evolvePDF(x2, mu, f2) ;

	double integrand = 0. ;
	for (int i=0; i < flav.size(); i++){
		// First we add contributions which don't require z integration
		// Contribution from Delta function
		integrand += 2* Iqfromq(1.,Lp).Delta * 
					f1[6+ flav[i] ] * f2[6 -flav[i] ] /tau/tau ; 
		// Contribution from plus function which is z independent 
		integrand += (PlusRem(x1) + PlusRem(x2)) * (-Lp/2) * 
						f1[6+ flav[i] ]*f2[6 -flav[i] ]/tau/tau ; 
	}

	// Contribution from remaining z dependent integration 
	double result, error ;
	zPara zpara;
	zpara.mu = mu; zpara.tau = tau; zpara.pTveto = pTveto;
	zpara.Y = Y; zpara.Nf = para.Nf ; zpara.flavor = para.flavor ;

	zpara.flag = 1; 
	NIntegrate(0., 1., &zpara, integrand_z, result, error, 0);
	integrand += result ;
	zpara.flag = 2; 
	NIntegrate(0., 1., &zpara, integrand_z, result, error, 0);
	integrand += result ;	

	return Jac*integrand ;
}

