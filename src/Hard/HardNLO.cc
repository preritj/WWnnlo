#include <iostream> 
#include <cstdlib> 
#include <vector> 
#include <math.h>
#include <const.h>
#include <input.h>
#include <PDF.h>
#include <Hard.h>
#include <misc.h>

using namespace std;
using namespace constants;

//----------------------------------------
//NLO hard function (process independent part)
double Hard::H1mu(double M, double mu_){
	set_mu(mu_); return H1mu(M) ;
}

double Hard::H1mu(double M){
	double L = 2.*log(mu/M) ;
	return  -CF*pi*pi/6 -CF* L*L - 3*CF*L ;  
}

//----------------------------------------
// NLO Hard function (process dependent part)
vector<flav> Hard::nlo(double M, double mu_){
	set_mu(mu_); return nlo(M); 
}

vector<flav> Hard::nlo(double M){
	vector<flav> H;
	flav H0 = lo(M) ;
	H.clear(); H.push_back(H0);
	flav result, error;
	flav H1;
	double H1mu_ = H1mu(M) ;
	hPara para(proc); 
	para.order= 1;
	para.M = M;

	//up-type 
	para.flavor = 'u';
	NIntegrate(0., 1., &para, integrand_cosT, result.u, error.u, 0);
	//up-type 
	para.flavor = 'd';
	NIntegrate(0., 1., &para, integrand_cosT, result.d, error.d, 0);

	H1.u = H1mu_*H0.u + result.u;
	H1.d = H1mu_*H0.d + result.d;
	H.push_back(H1);
	return H; 
}
