// z Integration
//-------------------------------------------------------------
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

double integrand_z(double x, void* para_ptr){
	double Jac = 1.0 ;  // Take care of Jacobian in the integrals
	zPara para = *(zPara*) para_ptr ;
	double tau = para.tau; double mu = para.mu; double Y=para.Y; 
    double pTveto = para.pTveto ; double Lp = 2*log(mu/pTveto) ;
	double x1 = para.tau*exp(Y); double x2 = para.tau*exp(-Y) ; 

	vector<int> flav; flav.clear();
	if (para.flavor == 'u') {flav.push_back(2); flav.push_back(4);} 
    else if (para.flavor == 'd') {
		flav.push_back(1); flav.push_back(3);
		if (para.Nf == 5) flav.push_back(5); }
	else {cout << "Wrong flavor type. Aborting..." << endl; abort();}
	
	double z_min = 0.;
	double z_max =  1. ;
	double z ;
	double f1[13], f2[13], f1p[13] ;
	if (para.flag == 1) {
		z_min = x1 ;
		z = z_min + x*(z_max - z_min) ;
		evolvePDF(x1/z, mu, f1) ; 
		evolvePDF(x1, mu, f1p) ; 
		evolvePDF(x2, mu, f2) ;
	}	
	else if (para.flag == 2) {
		z_min = x2 ;
		z = z_min + x*(z_max - z_min) ;
		evolvePDF(x2/z, mu, f1) ; 
		evolvePDF(x2, mu, f1p) ;
		evolvePDF(x1, mu, f2) ;
	}
	else {cout << "Invalid flag in NLO beam function. Aborting.. "<< endl;
			abort(); } 
	Jac = Jac * (z_max - z_min) ;

	double integrand = 0. ;
	for (int i=0; i < flav.size(); i++){
		int fl1 = 6 - (2*para.flag-3)*flav[i] ;
		int fl2 = 6 + (2*para.flag-3)*flav[i] ;
		// contribution from plus distribution
		integrand += Iqfromq(z,Lp).Plus * 
					(f1[fl1]-f1p[fl1])* f2[fl2] /tau/tau ; 
		integrand += Iqfromg(z,Lp).Plus * 
					(f1[6]-f1p[6])* f2[fl2] /tau/tau ; 
		//remaning contribution 
		integrand += Iqfromq(z,Lp).None * f1[fl1]* f2[fl2] /tau/tau ; 
		integrand += Iqfromg(z,Lp).None * f1[6] * f2[fl2] /tau/tau ; 
	}
	return integrand*Jac ;
}
