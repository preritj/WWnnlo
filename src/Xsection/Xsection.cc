#include <Xsection.h>
#include <PDF.h>

using namespace constants;
using namespace beamPara;

double Xsection::factor(double M){
	double s = M*M;
	vector<double> Mext = proc.ext_masses();
	double M1 = Mext[0];
	double M2 = Mext[1];
	double dr = (M1*M1 - M2*M2)/s ;
	double beta = sqrt(1.0 + dr*dr - 2/s*(M1*M1+M2*M2) ) ;
	return beta/16/pi/M/ECM/ECM*pbarn ;
}

void Xsection::set_mu(double mu_){
	mu = mu_; as = alpha_s(mu) ;
	hard.set_mu(mu_) ;
	beam.set_mu(mu_) ;
}

flav Xsection::lo(double M, double mu_){
	set_mu(mu_); return lo(M);
}

flav Xsection::lo(double M){
	double fac = factor(M);
	flav H0 = hard.lo(M) ;
	flav B0 = beam.lo(M) ;
	return 2*fac*H0*B0 ;
}

vector<flav> Xsection::nlo(double M, double mu_){
	set_mu(mu_); return nlo(M);
}

vector<flav> Xsection::nlo(double M){
	vector<flav> xs; xs.clear();
	double fac = factor(M);
	vector<flav> H = hard.nlo(M);
	flav B0 = beam.lo(M);
	flav B1 = beam.nlo(M);

	//LO xsection :
	xs.push_back(2*fac*H[0]*B0);

	// NLO xsection :
	flav Hnlo = H[0] + as/2/pi*H[1];
	xs.push_back(2*fac*(Hnlo*B0 + as/4/pi*H[0]*B1) );	
	return xs;
}
