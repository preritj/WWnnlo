

#include <iostream> 
#include <fstream> 
#include <complex>
#include <stdlib.h> 
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <vector>
#include <gsl/gsl_sf_zeta.h>
#include <gsl/gsl_sf_psi.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_dilog.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_vegas.h>
using namespace std;

//==========================================================================
// timer 
class timer
{
	private :
		clock_t start_time, stop_time, duration_time ;
	public : 
		void start() {start_time = clock() ; }
		void stop() {stop_time = clock() ; }
		double duration() {
			duration_time = stop_time - start_time ; 
			return double(duration_time)/((double)CLOCKS_PER_SEC) ;
		}
};

//==========================================================================
//Linking to PDFs (Fortran code for lhapdf) 
extern "C" 
{ 
//extern void PDFini_(); 
 void pdfini_(int&); 
 void evolvepdf_(double&, double&, double*  );
 double alphaspdf_(double&) ;
} 

// Returns PDF in f
void evolvePDF(double x, double Q, double* f)
{
	evolvepdf_(x , Q , f) ;
}

//running QCD coupling
double alpha_s (double mu2)
{
	double mu = sqrt(mu2) ;
	return alphaspdf_(mu) ;
}

//==========================================================================
//constants : SM parameters
const double pi = 3.14159265359 ;
const double zeta3 = gsl_sf_zeta_int(3), zeta5 = gsl_sf_zeta_int(5) ;
const double CA = 3.0, CF = 4.0/3, TF = 0.5, NC = 3; 
const double Euler = - gsl_sf_psi(1.0) ;
const double alpha = 1.0/130 ; 
const double pbarn = 3.894e8 ;
//const double MW = 80.376, MW2 = MW*MW ;
const double MW = 80.398, MW2 = MW*MW ;
const double MZ = 91.1876, MZ2 = MZ*MZ ;
//const double SW2 = 0.2312, SW = sqrt(SW2) ;
const double SW2 = 1.-MW2/MZ2, SW = sqrt(SW2) ;
const double GF = 1.16639e-5 ;
const int Nf = 5 ;

// constants that appear in beta function
const double B0 = (11*CA/3 - 4*TF*Nf/3) ;
const double B1 = (34*pow(CA,2)/3 - 20*CA*TF*Nf/3 - 4*CF*TF*Nf) ;
const double B2 = 2857*pow(CA,3)/54 
				+ (2*pow(CF,2) - 205*CF*CA/9 - 1415*pow(CA,2)/27)*
				TF*Nf + (44*CF/9 + 158*CA/27)*pow(TF*Nf,2) ;
// constants that appear in cusp anomalous dimension
const double G0 = (4*CF) ;
const double G1 = 4*CF*((67.0/9 - pow(pi,2)/3)*CA - 20.0/9*TF*Nf) ;
const double G2 = 4*CF*(pow(CA,2)*(245.0/6 -134*pow(pi,2)/27 +
					11*pow(pi,4)/45 + 22*zeta3/3) + CA*TF*Nf*
					(-418.0/27 + 40*pow(pi,2)/27 -56*zeta3/3) + 
					CF*TF*Nf*(-55.0/3+16*zeta3)-16*pow(TF*Nf,2)/27);
// constants that appear in quark anomalous dimension
const double gV0 = (-6*CF) ;
const double gV1 = pow(CF,2)*(-3.0 +4*pow(pi,2) -48*zeta3) + 
				+ CF*CA*(-961.0/27 -11*pow(pi,2)/3+52*zeta3)
				+ CF*TF*Nf*(260.0/27+4*pow(pi,2)/3);
const double gV2 = pow(CF,3)*(-29.0 - 6*pow(pi,2) - 16*pow(pi,4)/5 
					- 136*zeta3 + 32*pow(pi,2)/3*zeta3 + 480*zeta5 ) 
			+pow(CF,2)*CA*(-151.0/2 + 410*pow(pi,2)/9+494*pow(pi,4)/135
				 - 1688*zeta3/3 - 16*pow(pi,2)*zeta3/3 - 240*zeta5)
			+CF*pow(CA,2)*(-139345.0/1458 -7163*pow(pi,2)/243 -272*zeta5
			   - 83*pow(pi,4)/45 + 7052*zeta3/9 - 88*pow(pi,2)/9*zeta3)
			+pow(CF,2)*TF*Nf*(5906.0/27 -52*pow(pi,2)/9 -56*pow(pi,4)/27
				     + 1024*zeta3/9) 
			+CF*CA*TF*Nf*( -34636.0/729 + 5188*pow(pi,2)/243 
				      + 44*pow(pi,4)/45 - 3856*zeta3/27)
			+CF*pow(TF*Nf,2)*(19336.0/729-80*pow(pi,2)/27-64*zeta3/27) ;
// constants that appear in jet-veto functions
const double cLA = 131./72 - pi*pi/6 - 11./6*log(2.) ;
const double c0A = -805./216 + 11*pi*pi/72 + 35./18*log(2.)
					+11./6*pow(log(2.),2) + zeta3/2 ;
const double c2A = 1429./172800 + pi*pi/48 + 13./180*log(2.) ;
const double c4A = -0.0225794 ;
const double cLf = -23./36 + 2./3*log(2.) ;
const double c0f = 157./108 - pi*pi/18 - 8./9*log(2.) 
					- 2./3*pow(log(2.),2) ;
const double c2f = 3071./86400 - 7./360*log(2.) ;
const double c4f = -0.000442544 ;

//Non=perturbative scale
//const double LambdaNP = 0.0 ;

//==========================================================================
// Define complex number type : 
typedef complex<double> dcomplex;
// Overload power  function for complex type
dcomplex pow (dcomplex base, int exponent)
{
	dcomplex product = 1.0 ;
	for (int i = 1 ; i <= exponent ; i++) 
	product *= base;
	return product ;
}

//==========================================================================
// Define Kronecker Delta function
int Delta(int i, int j)
{
	if (i == j) return 1 ;
	else return 0 ;
}

//==========================================================================
// Numerical integration using gsl library 
// lower_limit = lower limit of integration
// upper_limit = upper limit of integration
// Ndim = # of integration variable
// params_ptr = pointer to parameters in the integration
// integrand( .. ) = function to be integrated 
// result = the result from the integration
// error = the estimated error from the integration
// Ncalls = Number of calls made for MC
// verbose = 1 : output results, 0 : no output
// VEGAS (MC) method if used when Ndim is explicitly specified
void NIntegrate(double* lower_limit, double* upper_limit, 
		int Ndim, void *params_ptr, double (*integrand)
		(double *, size_t, void*), double& result, 
		double& error, int Ncalls, int verbose )
{
	//gsl_integration_workspace *work_ptr 
	//	= gsl_integration_workspace_alloc (10000);

	double abs_error = 1.0e-3;	/* to avoid round-off problems */
	double rel_error = 1.0e-3;	/* the result will usually be much better */

	gsl_monte_function My_function;

	My_function.f = integrand;
	My_function.dim = Ndim;
	My_function.params = params_ptr;

	const gsl_rng_type *T;
    gsl_rng *r;
    gsl_rng_env_setup ();
    T = gsl_rng_default;
    r = gsl_rng_alloc (T);
    size_t calls = Ncalls;
	gsl_monte_vegas_state *s = gsl_monte_vegas_alloc (Ndim);
	
	gsl_monte_vegas_integrate (&My_function, lower_limit, 
	 upper_limit, Ndim, calls, r, s, &result, &error) ;
	 if (verbose == 1)
	 printf ("result = % .6f sigma = % .6f "
                "chisq/dof = %.1f\n", result, error, 
                gsl_monte_vegas_chisq (s));
    double res0, res_rel ;
    res0 = result ;
    res_rel = 1. ;
	
	do
      {
        gsl_monte_vegas_integrate (&My_function, lower_limit, 
	 	upper_limit, Ndim, calls, r, s, &result, &error);
	 	if (verbose == 1)
        printf ("result = % .6f sigma = % .6f "
                "chisq/dof = %.1f\n", result, error, 
                gsl_monte_vegas_chisq (s));
        res_rel = abs(result-res0)/result ;
        res0 = result ;
      }
    //while (fabs (error/result) > 0.001);
    while (fabs (gsl_monte_vegas_chisq (s) - 1.0) > 0.5
    			|| res_rel > 0.005);
	
	gsl_monte_vegas_free (s);
	gsl_rng_free (r);
} 

// QNG (quadrature) method is used if Ndim is not specified
// and it is assumed that Ndim = 1
void NIntegrate(double lower_limit, double upper_limit, 
		 void *params_ptr, double (*integrand)
		(double, void*), double& result, 
		double& error, int verbose )
{
	gsl_error_handler_t * old_handler=gsl_set_error_handler_off();
	gsl_integration_workspace *w 
		= gsl_integration_workspace_alloc (50000);

	double abs_error = 5.0e-3;	/* to avoid round-off problems */
	double rel_error = 1.0e-3;	/* the result will usually be much better */
	size_t Neval ;

	gsl_function My_function;

	My_function.function = integrand;
	My_function.params = params_ptr;
	
	int status=1;
    while(status) {
        status= gsl_integration_qags (&My_function, lower_limit, 
	 	upper_limit, abs_error, rel_error, 10000, w, &result, &error);
	 	if (status) rel_error*=2.0 ;
	 	if (verbose == 1) {
	 		printf ("result = % .6f error = % .6f \n", result, error);
	 		if (status) cout << "Increasing tolerance " << endl ;
	 		} ;
	 	}
    gsl_integration_workspace_free (w);
    gsl_set_error_handler(old_handler);
} 

//==========================================================================
// Define alpha_s for negative scale-squared
dcomplex alphas (double mu2){
	if (mu2 > 0.) return alpha_s(mu2) ;
	else {
		dcomplex a (0. , B0*alpha_s(-mu2)/4.) ;
		return alpha_s(-mu2)/(1. - a + B1/B0 * 
		     alpha_s(-mu2)/4./pi*log(1.- a)) ;
	} 
}

//==========================================================================
// Functions required for resummation
dcomplex agV (double nu2, double mu2, int Nloop)
{
	if (Nloop == 0) {
		dcomplex nu2eps (nu2,  - 1e-8) ;
		return -alphas(nu2)*gV0/8./pi*log(mu2/nu2eps) ;
	}
	else {
	dcomplex agV = 
		(gV0*log(alphas(mu2)/alphas(nu2)))/(2.*B0) ;
	if (Nloop > 1) {
		agV = agV + 
		((-(B1*gV0) + B0*gV1)*(alphas(mu2) - alphas(nu2)))/
   		(8.*(B0*B0)*pi) ; }
    if (Nloop > 2 ){
    	agV = agV + 
   		((B1*B1*gV0 - B0*B1*gV1 + B0*(-(B2*gV0) + B0*gV2))*
     	(pow(alphas(mu2),2) - pow(alphas(nu2),2)))/
   		(64.*pow(B0,3)*(pi*pi)) ; }
   	return agV ;
   	}
}

dcomplex aG (double nu2, double mu2, int Nloop)
{
	if (Nloop == 0) {
		dcomplex nu2eps (nu2,  - 1e-8) ;
		return -alphas(nu2)*G0/8./pi*log(mu2/nu2eps) ;
	}
	else {
	dcomplex aG = 
		(G0*log(alphas(mu2)/alphas(nu2)))/(2.*B0) ;
	if (Nloop > 1) {
		aG = aG + 
		((-(B1*G0) + B0*G1)*(alphas(mu2) - alphas(nu2)))/
   		(8.*(B0*B0)*pi) ; }
    if (Nloop > 2 ){
    	aG = aG + 
   		((B1*B1*G0 - B0*B1*G1 + B0*(-(B2*G0) + B0*G2))*
     	(pow(alphas(mu2),2) - pow(alphas(nu2),2)))/
   		(64.*pow(B0,3)*(pi*pi)) ; }
   	return aG ;
   	}
}

dcomplex SG (double nu2, double mu2, int Nloop)
{
	dcomplex r = alphas(mu2)/alphas(nu2) ;
	if (Nloop == 0) {
		dcomplex nu2eps (nu2,  -1e-8) ;
		return -G0*alphas(nu2)/32./pi*pow(log(mu2/nu2eps),2) ;
	}
	else {
	dcomplex SG = 
		-((G0*pi*(1. - r + r*log(r)))/(B0*B0*r*alphas(nu2))) ;
	if (Nloop > 1) {
		SG = SG + 
		(2.*(B1*G0 - B0*G1)*(-1. + r) + 
     	log(r)*(-2*B1*G0 + 2*B0*G1 + B1*G0*log(r)))/
  		 (8.*pow(B0,3)) ; }
    if (Nloop > 2 ){
    	SG = SG + 
   		-(alphas(nu2)*((B0*(B2*G0*(-1. - r) + 
              B0*G2*(-1. + r)) - B0*B1*G1*(-3. + r) + 
           B1*B1*G0*(-1. + r))*(-1. + r) + 
       	 2.*(B0*B2*G0 + B1*B1*G0*(-1. + r) - 
           B0*B1*G1*r)*log(r)))/(32.*pow(B0,4)*pi) ; }
   	return SG ;	
   	}
}

dcomplex WilsonC_evolve (double mu2_i, double mu2_f, 
						double M2, int Nloop)
{
	dcomplex M2eps (-M2, - 1e-8) ;
	if (Nloop == 0) {
		return (2.*SG(mu2_i, mu2_f, 0)-aG(mu2_i, mu2_f, 0)* 
		log(M2eps/mu2_i) - agV(mu2_i, mu2_f, 0)) ;
	}
	else 
	return exp(2.*SG(mu2_i, mu2_f, Nloop+1)-aG(mu2_i, mu2_f, Nloop)* 
		log(M2eps/mu2_i) - agV(mu2_i, mu2_f, Nloop) ) ;
}


double Fqq(double pT_veto, double mu, double R, int Nloop)
{
	double LT = 2*log(mu/pT_veto) ; 
	double as = alpha_s(mu*mu)/4/pi ;
	double Fqq = as*G0*LT ;

	if (Nloop > 1) {
		double fR = CA*(cLA * log(R) + c0A + c2A*R*R + c4A*pow(R,4)) 
					+ CF*(-pi*pi*R*R/12 + pow(R,4)/16)
					+ TF*Nf*(cLf*log(R) +c0f +c2f*R*R + c4f*pow(R,4)) ;
		double d2 = CF*((808./27 - 28*zeta3)*CA - 224./27*TF*Nf);
		d2 = d2 - 32*CF*fR ;
//***********************************BULLSHIT******************************
//		d2 = 0. ;
//*************************************************************************
		
		Fqq = Fqq + as*as*(G0*B0*LT*LT/2 + G1*LT + d2) ;
	}

	return Fqq ;
}

double hF(double pT_veto, double mu)
{
	double LT = 2*log(mu/pT_veto) ; 
	double as = alpha_s(mu*mu)/4/pi ;
	return as*(-G0*LT*LT/4 - gV0/2*LT) ;
}

//==========================================================================
// Functions for LO and NLO amplitudes / Wilson coefficients
double C_SS_u (double S)
{
	double C_SS = 
	(9*(S*S) - 24*MZ2*S*SW2 + 32*(MZ2*MZ2)*(SW2*SW2))/
    (3.*pow(MZ2 - S,2)) ;
    return C_SS ;
}

double C_SS_d (double S)
{
	double C_SS = 
	(9*(S*S) - 12*MZ2*S*SW2 + 8*(MZ2*MZ2)*(SW2*SW2))/
	(3.*pow(MZ2 - S,2));
    return C_SS ;
}

double C_ST_u (double S)
{
	double C_ST = 
	(3*S - 4*MZ2*SW2)/(MZ2 - S) ;
    return C_ST ;
}

double C_ST_d (double S)
{
	double C_ST = 
	(3*S - 2*MZ2*SW2)/(MZ2 - S) ;
    return C_ST ;
}

double F0_SS (double S, double T, double U)
{
	double F0 = 
	(-4*pow(MW2,4) + 4*MW2*S*(S*S + S*T + T*T) - 
     MW2*MW2*(17*(S*S) + 20*S*T + 12*(T*T)) + 
     4*pow(MW2,3)*(5*T - U) + S*S*T*U)/(MW2*MW2*(S*S)) ;
     return F0 ;
}

double F0_TT (double S, double T, double U)
{
	double F0 = 
	(3.*(-4*pow(MW2,4) + 4*MW2*S*(T*T) - MW2*MW2*T*(T - 4*U) + 
       pow(T,3)*U))/(MW2*MW2*(T*T)) ;
     return F0 ;
}

double F0_ST (double S, double T, double U)
{
	double F0 = 
	(4*pow(MW2,4) + pow(MW2,3)*(8*S - 6*T) - 5*(MW2*MW2)*S*T + 
     2*MW2*T*(2*(S*S) + S*T + T*T) + S*(T*T)*U)/(MW2*MW2*S*T) ;
     return F0 ;
}


dcomplex F1_SS (double S, double T, double U, double mu2, int pole)
{
	double theta = 1. ;
	if (mu2 < 0.) {mu2 = - mu2 ; theta = 0. ;} 
	double finite, pole_1, pole_2 ;
	if (pole == 0) {
		finite = 1. ; pole_1 = 0. ; pole_2 = 0; 
	}
	else if (pole == -1)  {
		finite = 0. ; pole_1 = 1. ; pole_2 = 0; 
	}
	else if (pole == -2)  {
		finite = 0. ; pole_1 = 0. ; pole_2 = 1.; 
	}
	//double Ls = 0.; 
	dcomplex Ls ( log(S/mu2) , - pi*theta ) ;
	dcomplex C1 = 1/S*(pole_2 - pole_1*Ls + 
	    finite* Ls*Ls/2.) ;
		
	dcomplex F1 = (finite*(8./3 - Ls) + pole_1 + 2./3*C1*S)*
	(4*pow(MW2,4) - 4*MW2*S*(S*S + S*T + T*T) + 
     MW2*MW2*(17*(S*S) + 20*S*T + 12*(T*T)) - S*S*T*U 
     + 4*pow(MW2,3)*(-5*T + U)
     )/(MW2*MW2*(S*S)) ;
     

	return F1 ;
}

dcomplex F1_TT (double S, double T, double U, double mu2, int pole)
{
	double theta = 1. ;
	if (mu2 < 0.) {mu2 = - mu2 ; theta = 0. ;} 
	double finite, pole_1, pole_2 ;
	if (pole == 0) {
		finite = 1. ; pole_1 = 0. ; pole_2 = 0; 
	}
	else if (pole == -1)  {
		finite = 0. ; pole_1 = 1. ; pole_2 = 0; 
	}
	else if (pole == -2)  {
		finite = 0. ; pole_1 = 0. ; pole_2 = 1.; 
	}
	
	double  betaW = sqrt(1.0 -4*MW2/S) ;
	//double Ls = 0. ; 
	dcomplex Ls ( log(S/mu2) , - pi*theta ) ;
	//double Lt = log(-T/S) ; 
	dcomplex Lt ( log(-T/mu2) , pi*(1.-theta) ) ;
	//double Lm = log(MW2/S) ; 
	dcomplex Lm ( log(MW2/mu2) , -pi*theta ) ;
	
	dcomplex C1 = 1/S*(pole_2 - pole_1*Ls + 
	    finite* Ls*Ls/2. ) ;
	dcomplex C2 = (pole_1*(Lt-Lm) + finite*(
	     Lm*Lm/2. - Lt*Lt/2.))/(MW2-T) ;
	dcomplex C3 = finite/betaW/S*(4*gsl_sf_dilog((betaW-1)/(betaW+1))
				+pi*pi/3  + pow(log((1-betaW)/(1+betaW)),2) ) ;
	dcomplex D1 = 1/S/T*(pole_2 - pole_1*(2.*Lt + Ls - 2.*Lm) + 
	      finite*(2*pi*pi/3+ Ls*Ls/2. + 2.*Ls*(Lt-Lm)  + 
				2.*pow(Lt-Lm,2) -4*gsl_sf_dilog(T/MW2)
				-4*log(1-T/MW2)*(Lt-Lm) )) ;

     
	 dcomplex F1 = finite*(
	40. + (36*(MW2*MW2))/(T*T) - (70*MW2)/T + (36*S)/T - 
   (2*pow(MW2 + T,2))/(betaW*betaW*MW2*T) + 
   (-32*S + 2*T)/MW2 - (8*T*U)/(MW2*MW2)) ;
    
    F1 = F1 - (finite*Ls-pole_1)*(
    7. - (5*MW2)/T - (3*pow(MW2 + T,2))/(pow(betaW,4)*MW2*T) + 
   (4*(-4*(MW2*MW2) - 4*MW2*T + T*T))/(betaW*betaW*MW2*T) - 
   (3*T*U)/(MW2*MW2) + (-11*S + U)/MW2) ; 
   
   F1 = F1 - (finite*Lt-pole_1)*(
   4.*(1 - (3*(MW2*MW2))/(T*T) + (S*(3*MW2 + T))/(T*(-MW2 + T))) ) ;
     
    F1 = F1 - (finite*Lm-pole_1)*(
    2. + (16*S)/(MW2 - T) + (24*(MW2*MW2))/(T*T) + (-19*MW2 
    + 24*S)/T + T/MW2 + (3*pow(MW2 + T,2))/(pow(betaW,4)*MW2*T) + 
   (4*(4*(MW2*MW2) + 4*MW2*T - T*T))/(betaW*betaW*MW2*T)) ; 
   
   F1 = F1 + C2*(
   (8*(MW2 - T)*(2*(MW2*MW2) - T*U))/(T*T) ) ;
   
   F1 = F1 + C1*(
   2*S*(1. - (2*U)/T - (4*MW2*S + T*U)/(MW2*MW2))) ;
   
   F1 = F1 + C3*(
   2*MW2 - 4*S - (3*(MW2*MW2) + 4*(S*S))/T + T - 
   (3*pow(MW2 + T,2))/(pow(betaW,4)*T) + 
   (2*(-9*(MW2*MW2) - 10*MW2*T + T*T))/(betaW*betaW*T)) ;
   
   F1 = F1 + D1*(
   (2*S*(4*MW2*S - S*S + T*T + U*U))/T ) ; 

	return F1 ;
}

dcomplex F1_ST (double S, double T, double U, double mu2, int pole)
{
	double theta = 1. ;
	if (mu2 < 0.) {mu2 = - mu2 ; theta = 0. ;} 
	double finite, pole_1, pole_2 ;
	if (pole == 0) {
		finite = 1. ; pole_1 = 0. ; pole_2 = 0; 
	}
	else if (pole == -1)  {
		finite = 0. ; pole_1 = 1. ; pole_2 = 0; 
	}
	else if (pole == -2)  {
		finite = 0. ; pole_1 = 0. ; pole_2 = 1.; 
	}
	
	double  betaW = sqrt(1.0 -4*MW2/S) ;
	//double Ls = 0. ; 
	dcomplex Ls ( log(S/mu2) , - pi*theta ) ;
	//double Lt = log(-T/S) ; 
	dcomplex Lt ( log(-T/mu2) , pi*(1.-theta) ) ;
	//double Lm = log(MW2/S) ; 
	dcomplex Lm ( log(MW2/mu2) , -pi*theta ) ;
	
	dcomplex C1 = 1/S*(pole_2 - pole_1*Ls + 
	    finite* Ls*Ls/2. ) ;
	dcomplex C2 = (pole_1*(Lt-Lm) + finite*(
	     Lm*Lm/2. - Lt*Lt/2.))/(MW2-T) ;
	dcomplex C3 = finite/betaW/S*(4*gsl_sf_dilog((betaW-1)/(betaW+1))
				+pi*pi/3  + pow(log((1-betaW)/(1+betaW)),2) ) ;
	dcomplex D1 = 1/S/T*(pole_2 - pole_1*(2.*Lt + Ls - 2.*Lm) + 
	      finite*(2*pi*pi/3+ Ls*Ls/2. + 2.*Ls*(Lt-Lm)  + 
				2.*pow(Lt-Lm,2) -4*gsl_sf_dilog(T/MW2)
				-4*log(1-T/MW2)*(Lt-Lm) )) ;
				
	dcomplex F1 = finite*(
	(-4.*(17*pow(MW2,5) + pow(MW2,4)*(34*S - 43*T) + 
       pow(MW2,3)*T*(-55*S + 27*T) - 
       4*MW2*(T*T)*(5*(S*S) + 3*S*T + 2*(T*T)) + 
       MW2*MW2*T*(18*(S*S) + 37*S*T + 7*(T*T)) - 4*S*pow(T,3)*U)
     )/(3.*(MW2*MW2)*S*(MW2 - T)*T)) ; 
   
   F1 = F1 - (finite*Ls-pole_1) * (
   (-4.*(MW2 + T))/(3.*(betaW*betaW)*MW2) - 
   (2.*(6*pow(MW2,4) + 6*pow(MW2,3)*(2*S - 3*T) + 
        2*MW2*T*(6*(S*S) + 2*S*T + 3*(T*T)) + 3*S*(T*T)*U + 
        MW2*MW2*T*(-5*S + 6*U)))/(3.*(MW2*MW2)*S*T) ) ;
        
    F1 = F1 - (finite*Lt-pole_1) * (
    (4.*(3*pow(MW2,4) + pow(MW2,3)*(6*S - 8*T) + 
       MW2*MW2*T*(-7*S + 6*T) + 
       T*(2*MW2 - U)*(-(T*T) + 2*S*(S + U))))/
   (3.*S*pow(MW2 - T,2)*T) ) ;
   
   F1 = F1 - (finite*Lm-pole_1) * (
   (4.*(MW2 + T))/(3.*(betaW*betaW)*MW2) + 
   (4.*(5. + (8*MW2)/S - (2*MW2*S)/pow(MW2 - T,2) - 
        (6*MW2*(MW2 + 2*S))/(S*T) - ((2*MW2 + S)*T)/(MW2*S) + 
        (2*(4*MW2 + S))/(-MW2 + T)))/3. ) ;

	F1 = F1 + C2 * (
	(-8.*(MW2 - T)*(2*MW2*(MW2 + 2*S) - (2*MW2 + S)*T + T*T))/
	(3.*S*T) ) ;
   
    F1 = F1 + C1 * (
    (-4.*(2*pow(MW2,4) + pow(MW2,3)*(4*S - 6*T) + 
       2*MW2*T*(2*(S*S) + S*T + T*T) + S*(T*T)*U + 
       MW2*MW2*T*(-3*S + U)))/(3.*(MW2*MW2)*T) ) ;
       
    F1 = F1 + C3 * (
    (-4.*(MW2 + T))/(3.*(betaW*betaW)) - 
   (2.*(8*(S*S) + 3*S*T + T*T + 5*S*U - U*U))/(3.*S) ) ;
   
   F1 = F1 + D1 * (
   (-2.*(4*MW2*S + T*T + (S + U)*(3*S + U)))/3. ) ;
   
   return F1 ;

}

template<typename Type_parton>
struct partons
{
	Type_parton u, d, g, q, uu, dd, ug, dg ;
} ;

//==========================================================================
class xsection
{
  private : 
	int Nf, Nloop ;  // Nf = Active flavors, Nloop = QCD Order for resum
	int part; // relevant for fixed order : 
	// part 0 and part 1 should be added separately
	double ECM; // ECM = Center of mass energy
	double M, M2, betaW, Y, cosT ; // Phase space variables
	double S, T, U ; // Mandelstam variables
	double mu_f2, mu_r2, mu_RG2, mu_h2 ; // scales
	double pT_veto, R ; // jet variables
	bool Matching ; // flag if set to true expands NLL results to fixed order
	partons<double> ff0 ; // LO Luminosity and squared amplitude
	double scaleFactor ; // scale factor for dynamic scale
	double NuScaleFactor ; // scale factor for rapidity scale
	bool dynamic_scale ; // if true : dynamic scale used for hard matching
	bool is_b; // if true : b-quark included 
	bool is_diff; // if true : output diff xsection
//	double mu_h, mu_s, mu_f ; // different scales entering factorization 
//	double G0,G1,G2,G3,gV0,gV1,gV2,gP0,gP1,gP2; //anomalous dimensions
//	double B0, B1, B2, B3 ; // Beta functions at various orders
//	double s, M, Y ; // E_CM, lepton pair invariant mass, rapidity
//	double eta, eta_p ;
//	double soft (double, double, double) ;
	// Initialization routine. Comes by default with the constructor
//	void init(); 
//	friend double integrand_xs_diffM(double, void* ) ;
//	friend double xsection_diffM(double, double, xsection) ;
//	friend double integrand_xs(double*, size_t, void* ) ;
	friend double integrand_diffM(double, void* ) ;
	
  public : 
  	xsection() ; // default constructor
  	xsection(int, int) ;
  	// routines to set Flavor number, scales and loop order 
//  	void set_Nf(int) , set_Nloop(int), set_scales(double, double, double);
  	//parton luminosity functions for differential xsections
//  	dcomplex lum_diffM (double, double) ;
  	// some functions for hard coefficient made public 
//  	dcomplex hard(double, double), hard_evolve(double, double, double);	
  	// Some coefficients appearing in hard and soft functions evolution : 
//	double aG(double, double), agV(double, double), agP(double, double) ;
//	double S(double, double), U(double, double, double, double) ;
//	double WilsonC (double, double) ; // Final Wilson coefficient 
//	double xsection_diffM (double, double) ; //diff xsection in M2
	int get_Nloop() ;
	bool isMatching() ;
	void set_PhaseSpace(double, double , double) ; //M,Y, cosT
	void set_Mandelstam(double, double , double) ; //S,T,U
	void set_scales() , set_ECM(double) , set_Match(bool) ;
	void set_muf2(double) , set_mur2(double), set_R(double) ;
	void set_muh2(double),set_muRG2(double),set_scaleFactor(double, double);
	//void set_muh2(double),set_muRG2(double),set_scaleFactor(double);
	void set_Nloop(int), incl_b(bool), set_diff(bool);
	double get_ECM() ;
	partons<double> Lum() , Hard() ;
	partons<dcomplex> WilsonC(double) ;
	double factor(), Col_anomaly() ;
	void set_pTveto(double), evaluate(), set_dynamic(bool);
} ; 

//==========================================================================
//  Function for initialization of the xsection class
xsection::xsection()
{
	Nf = 5;
	Nloop = 0;
	ECM = 8000.0 ;
	pdfini_(Nloop) ;
	mu_f2 = 0. ; mu_r2 =0.;
	mu_h2 = 0. ; mu_RG2 = 0.;
	dynamic_scale = false ;
	//NuScaleFactor = 1.0 ;
	pT_veto = 5. ;
	R = 0.4 ;
	Matching = false ;
	is_b = false ;
	is_diff = false ;
}

//==========================================================================
xsection::xsection(int N, int N_pdf)
{
	Nf = 5;
	Nloop = N;
	ECM = 8000.0 ;
	pdfini_(N_pdf) ;
	mu_f2 = 0. ; mu_r2 =0.;
	mu_h2 = 0. ; mu_RG2 = 0.;
	dynamic_scale = false ;
	//NuScaleFactor = 1.0 ;
	pT_veto = 25. ;
	R = 0.4 ;
	Matching = false ;
	is_b = false ;
	is_diff = false ;
} 

//  Function to set/retrieve loop
void xsection::set_Nloop(int N)
{ Nloop = N; }

int xsection::get_Nloop()
{ return Nloop; }

//  Function that determines if b-quark contribution is included
void xsection::incl_b(bool is_bq)
{ is_b = is_bq ;}

// Function that outputs differential distributions
void xsection::set_diff(bool isdiff)
{ is_diff = isdiff ;}

//==========================================================================
// Function to set/retrieve flag for Matching :
// If argument is false : NLL resummation is performed
// If argument is true : NLL results are expanded to fixed order
void xsection::set_Match(bool isMatching) 
{
	if (isMatching == true) Matching = true ;
	else Matching = false ;
}

bool xsection::isMatching() 
{ return Matching ;}

//==========================================================================
//  Function to set/retrieve CoM energy or sqrt{S}
void xsection::set_ECM(double ECM_in)
{
	ECM = ECM_in ;
}

double xsection::get_ECM() {return ECM ;}

//==========================================================================
//  Functions to set up phase-space for the xsection
void xsection::set_PhaseSpace(double M_in,  double Y_in, double cosT_in)
{
	M = M_in ; Y = Y_in ; cosT = cosT_in ; 
	M2 = M*M ; betaW = sqrt(1.0 -4*MW2/M2) ;
	
	//define Mandelstam variables in CoM frame
	double S = M2 ;
	double T = MW2 - M2/2* (1. - betaW * cosT);
	double U = 2*MW2 - S - T ;
	//set the Mandelstam variables for the xsection :
	set_Mandelstam(S, T, U) ;
	
	//set the scales
	set_scales() ;
}

//==========================================================================
//  Functions to set up Mandelstam variables for the xsection
void xsection::set_Mandelstam(double S_in, double T_in, double U_in)
{
	S = S_in ; T = T_in ; U = U_in ; 
}

//==========================================================================
//  Functions to set scales
void xsection::set_dynamic(bool isDynamic)
{ dynamic_scale = isDynamic ; }

void xsection::set_muf2(double muf2)
{ mu_f2 = muf2 ; }

void xsection::set_mur2(double mur2)
{ mu_r2 = mur2 ; }

void xsection::set_scaleFactor(double factor, double nuFactor)
//void xsection::set_scaleFactor(double factor)
{ scaleFactor = factor ;
	NuScaleFactor = nuFactor ;
}

void xsection::set_muh2(double muh2)
{ mu_h2 = muh2 ;}

void xsection::set_muRG2(double muRG2)
{ mu_RG2 = muRG2 ; }

void xsection::set_pTveto(double pTveto)
{ pT_veto = pTveto ; 
set_muRG2(pow(pT_veto,2)) ; }

void xsection::set_scales()
{
	if (abs(mu_f2) < 1.) set_muf2(MW2) ;
	if (abs(mu_r2) < 1.) set_mur2(MW2) ;
	if (abs(mu_h2) < 1.) set_muh2(MW2) ;
	if (abs(mu_RG2) < 1.) set_muRG2(MW2) ;
}

void xsection::set_R(double R_in)
{ R = R_in ; }
//==========================================================================
//  Function to return parton luminosity
// mu is the scale, u is the parameter at which luminosity is evaluated
//partons_ID = 0 : uu, 1 : dd, 2 : ug, 3 : dg
//dis_type  is the type of distribution that is returned :

//==========================================================================
class luminosity
{
  private : 
	double mu, tau, pT_veto;
	bool is_b;
	int partons_ID;
	partons<bool> flag ;
	friend double integrand_beam(double, void* ) ;
	
  public : 
    partons<double> ff0 ;
	partons<double> parton_lum (int , double) ;
	partons<double> kernels(int , double ) ;
	
	void set_para (double mu_in, double tau_in, 
		double pTveto_in, bool is_bin)
	{mu = mu_in; tau = tau_in; pT_veto = pTveto_in;
	 is_b=is_bin; ff0 = parton_lum(0, 1.) ;}
	
	void select_u() { flag.u = true ; flag.d = false ; }
	void select_d() { flag.d = true ; flag.u = false ; }
} ;

struct para_lum
{
	double tau, mu ;
	int partons_ID ;
	bool is_b ;
} ;

//  Integrand of the parton luminosity
// The limits for x are taken from 0 to 1 and the appropriate Jacobian included
double integrand_lum(double x, void* para_ptr)
{
	double eps = 1e-8 ;
	double Jac = 1.0 ;  // Take care of Jacobian in the integrals
	para_lum para = *(para_lum*) para_ptr ;
	double tau = para.tau ;
	double mu = para.mu ;
	int partons_ID = para.partons_ID ;
	int n_d = 2; //# of down-type quarks
	if (para.is_b) n_d = 3;
	
	double Y_min = - log(1/tau)/2. ;
	double Y_max =  log(1/tau)/2. ;
	double Y = Y_min + x*(Y_max - Y_min) ;
	Jac = Jac * (Y_max - Y_min) ;
	double z = sqrt(tau)*exp(Y) ;
	
	double f1[13], f2[13] ;
	evolvePDF(z, mu, f1) ;
	evolvePDF(tau/z, mu, f2) ;
	
	double integrand = 0. ;
	if (partons_ID == 0) {
		for (int i = 0 ; i < 2 ; i++) {	
			integrand = integrand + 2*f1[4-2*i]*f2[8+2*i];
		}	
	}
	else if (partons_ID == 1) {
		for (int i = 0 ; i < n_d ; i++) {	
			integrand = integrand + 2*f1[7+2*i]*f2[5-2*i];
		}	
	}
	else if (partons_ID == 2) {
		for (int i = 0 ; i < 2 ; i++) {	
			integrand = integrand 
					+ 2*f1[6]*(f2[4-2*i]+f2[8+2*i]) ;
		}	
	}
	else if (partons_ID == 3) {
		for (int i = 0 ; i < n_d ; i++) {	
			integrand = integrand 
					+ 2*f1[6]*(f2[7+2*i]+f2[5-2*i]) ;
		}	
	}
			
	return Jac*integrand/tau ;
}


//==========================================================================
//  Function to return product of beam kernels
// mu is the scale, pT_veto is the jet-veto scale
// z is the parameter at which beam kernels are evaluated 
//dis_type  is the type of distribution that is returned :
// 0 : delta function, 1 : plus distribution, 2 : the usual case  
// flag : flag.u = true evaluates up-type beam product, similarly for d


partons<double> luminosity::parton_lum (int dis_type, double z)
{
	double result, error; 
	double Lum[4] = {0.,0.,0.,0.} ;
	para_lum para ;
	para.mu = mu ; 
	para.is_b = is_b ;
	if ( dis_type == 0 ) {z = 1.; para.tau = tau ;}
	else para.tau = tau/z;
	
	
	//printf("tau = % .3f , mu = % .3f \n" , tau, mu) ;
	//abort() ;
	
	for (int pID = 0 ; pID < 4 ; pID++) {
		para.partons_ID = pID ;
		NIntegrate(0., 1., &para, 
				integrand_lum, result, error, 0);
		//printf("tau = % .3f , mu = % .3f , lum = %.5f \n" , tau, mu, result ) ;
		Lum[pID] = Lum[pID] +  result ;
		// Plus distribution
		if ( dis_type == 1 ) {
				para.tau = tau ;
				NIntegrate(0., 1., &para, 
						integrand_lum, result, error, 0);
				Lum[pID] = Lum[pID] - result ;
				para.tau = tau/z ;
			}
		}
	
	partons<double> lum ;
	lum.uu = Lum[0] ; lum.dd = Lum[1] ; 
	lum.ug = Lum[2] ; lum.dg = Lum[3] ;
	
	return lum ;
}

partons<double> luminosity::kernels(int dist_type, double z)
{
	double eps = 1e-8 ;
	double mu2 = mu*mu ;
	
	double LT = 2*log(mu/pT_veto) ; 
	double as = alpha_s(mu2)/4/pi ;
	partons<double> II;
	double Rqq, Rqg, Pqq, Pqg ;
	
 	if (dist_type == 0) {
 	    Pqq = 4*CF*(tau*(1+tau/2.) + 2.*log(1-tau)) ;
 		Pqg = 0. ;	
 		Rqq = - CF*pi*pi/6 ; 
 		Rqg = 0. ;	
 		} 
 	else if (dist_type == 1) {
 		Pqq = 4*CF*(1+z*z)/(1-z+eps) ;
 		Pqg = 0. ;	
 		Rqq = 0. ; 
 		Rqg = 0. ;		
 		} 
	else {
		Pqq = 0. ;
 		Pqg = 4*TF*(z*z + pow(1.-z,2)) ;	
 		Rqq = 2*CF*(1-z) ; 
 		Rqg = 4*TF*z*(1-z) ;		
	}
	
	II.q = as*(- Pqq*LT/2 + Rqq) ; 
 	II.g = as*(- Pqg*LT/2 + Rqg) ;		
	return II ;
}

//  Integrand of the beam product
// The limits for x are taken from 0 to 1 and the appropriate Jacobian included
double integrand_beam(double x, void* para_ptr)
{
	double eps = 1e-8 ;
	luminosity lum = *(luminosity*) para_ptr ;
	double mu = lum.mu ;
	double mu2 = mu*mu ;
	double tau = lum.tau ;
	double Jac = 1.0 ;  // Take care of Jacobian in the integrals
	
	double Y_min = - log(1/tau)/2. ;
	double Y_max =  log(1/tau)/2. ;
	double Y = Y_min + x*(Y_max - Y_min) ;
	Jac = Jac * (Y_max - Y_min) ;
	double z = sqrt(tau)*exp(Y) ;
	
	double integrand = 0. ;
	partons<double> II = lum.kernels(2,z);
	partons<double> IIplus = lum.kernels(1,z);
	partons<double> ff0 = lum.ff0 ;
	partons<double> ff =  lum.parton_lum(2,z);
	partons<double> ffplus = ff ;
	ffplus.uu = ffplus.uu - z*ff0.uu ; 
	ffplus.ug = ffplus.ug - z*ff0.ug ;
	ffplus.dd = ffplus.dd - z*ff0.dd ; 
	ffplus.dg = ffplus.dg - z*ff0.dg ;
	
	if (lum.flag.u == true) 
	integrand = 2*II.q*ff.uu + II.g*ff.ug 
			+ 2*IIplus.q*ffplus.uu + IIplus.g*ffplus.ug ;
	if (lum.flag.d == true) 
	integrand = 2*II.q*ff.dd + II.g*ff.dg 
			+ 2*IIplus.q*ffplus.dd + IIplus.g*ffplus.dg;
	
	return Jac*integrand ;
}

//==========================================================================
//  Functions for giving total luminosity i.e.
//	Nloop = 0   : product of PDFs
//  Nloop = 1   : product of Beam functions (fixed order)
//	Nloop = 2/3 : product of Beam functions
partons<double> xsection::Lum()
{
	double mu, mu2 ;
	double tau = M*M/ECM/ECM ;
	partons<double> Lum ;
	
	if (Nloop == 0) mu2 = mu_f2 ;
	else mu2 = mu_RG2 ;
	mu = sqrt(mu2) ;
	
	//printf("M = % .3f , tau = % .3f , mu = % .3f \n" , M, tau, mu) ;
	luminosity lum ;
	lum.set_para (mu, tau, pT_veto, is_b) ;
	partons<double> ff0 = lum.ff0 ;
	//cout << "test" ; abort() ;
	
	if (Nloop == 0 || Nloop == 2 || Nloop == 4
	    || (Nloop == 1 && part == 0)) 
	{Lum.u = ff0.uu ; Lum.d = ff0.dd ;}


	else {
		partons<double> II = lum.kernels(0,1.);
		Lum.u = 2*II.q*ff0.uu + II.g*ff0.ug ;
		Lum.d = 2*II.q*ff0.dd + II.g*ff0.dg ;
		if (Nloop >= 2) {
			Lum.u = Lum.u + ff0.uu ;
			Lum.d = Lum.d + ff0.dd ;
		}
		
		double result, error ;
		lum.select_u() ;
		NIntegrate(0., 1., &lum, integrand_beam, result, error, 0);
		Lum.u = Lum.u + result ;
		
		lum.select_d() ;
		NIntegrate(0., 1., &lum, integrand_beam, result, error, 0);
		Lum.d = Lum.d + result ;
	}
	
	return Lum ;
}


class hard
{
  private : 
	int Nloop ;
	double M2, S, T, U, mu2 ;
	partons<double> AmpSq0 ;
	partons<bool> flag ;
	
  public :
    partons<dcomplex> WilsonC() ;
    partons<double> AmpSq();
    void set_Mandelstam(double) ;
    void AmpSq_init() ;
    
    void set_para (double mu2_in, double M2_in, int N_in)
	{mu2 = mu2_in; M2 = M2_in; Nloop = N_in;}
	
	void select_u() { flag.u = true ; flag.d = false ; }
	void select_d() { flag.d = true ; flag.u = false ; }
	
	bool IsParton_u() { return (flag.u == true) ; }
	bool IsParton_d() { return (flag.d == true) ; }
};

//==========================================================================
//  Functions for giving Wilson coefficients 
partons<dcomplex> hard::WilsonC()
{		
	partons<dcomplex> WilsonC;
	dcomplex as = alphas(mu2)/pi/2. ;
	dcomplex AmpSq ;
	
	dcomplex f_SS = F1_SS(S,T,U,mu2,0) ;
			//- pi*pi/12*F1_SS(S,T,U,mu2,-2) 
	dcomplex f_ST = F1_ST(S,T,U,mu2,0)  ;
			//- pi*pi/12*F1_ST(S,T,U,mu2,-2) 
	dcomplex f_TT = F1_TT(S,T,U,mu2,0) ;
			//- pi*pi/12*F1_TT(S,T,U,mu2,-2)  
			
	AmpSq = 2.*f_SS * C_SS_u(S) 
			+ 2.*f_ST * C_ST_u(S) 
			+ 2.*f_TT ;
	WilsonC.u = 1. + as*(AmpSq/AmpSq0.u + pi*pi*G0/48.) ;
	
	AmpSq = 2.*f_SS * C_SS_d(S) 
			+ 2.*f_ST* C_ST_d(S) 
			+ 2.*f_TT ;
	WilsonC.d = 1. + as*(AmpSq/AmpSq0.d + pi*pi*G0/48.) ;
	
	WilsonC.g = 1. ; 
	
	return WilsonC ;
}

//==========================================================================
//  Functions for amplitude squared
void hard::AmpSq_init()
{	
	partons<double> AmpSq;
	AmpSq.u = F0_SS(S,T,U) * C_SS_u(S) 
			+ 2*F0_ST(S,T,U) * C_ST_u(S) + F0_TT(S,T,U);
	AmpSq.d = F0_SS(S,T,U) * C_SS_d(S) 
			+ 2*F0_ST(S,T,U) * C_ST_d(S) + F0_TT(S,T,U) ;
		
	AmpSq0 = AmpSq;	
}

partons<double> hard::AmpSq()
{	
	partons<double> AmpSq;
	if (Nloop == 0 || Nloop == 2) {
		AmpSq.u = AmpSq0.u;
		AmpSq.d = AmpSq0.d ;
	}
	else if (Nloop == 1) {
		AmpSq.u = AmpSq0.u*real(2.*WilsonC().u -1.) ;
		AmpSq.d = AmpSq0.d*real(2.*WilsonC().d -1.) ;
	}
	else if (Nloop > 2) {
		AmpSq.u = AmpSq0.u*real(2.*WilsonC().u-1.) ;
		AmpSq.d = AmpSq0.d*real(2.*WilsonC().d-1.)  ;
	}
		
	return AmpSq;	
}

//==========================================================================
//  Functions to set up Mandelstam variables for the xsection
void hard::set_Mandelstam(double cosT)
{
	double betaW = sqrt(1.0 -4*MW2/M2) ;
	
	//define Mandelstam variables in CoM frame
	 S = M2 ;
	 T = MW2 - M2/2* (1. - betaW * cosT);
	 U = 2*MW2 - S - T ;
	 AmpSq_init() ;
}

double integrand_cosT(double x, void* para_ptr)
{
	hard H = *(hard*) para_ptr ;
	
	double Jac = 1.0 ;  // Take care of Jacobian in the integrals
	// Integration over cos(theta) 
	// theta is the angle between p3 and the x-axis in CoM frame
	double cosT_min = -1.0 ;
	double cosT_max = 1.0 ;
	double cosT = cosT_min + x*(cosT_max - cosT_min) ;
	Jac = Jac*(cosT_max - cosT_min) ; 
	
	H.set_Mandelstam(cosT) ;
	//printf("S = % .3f , T = % .3f , U = % .3f \n" , S,T,U) ;
	//cout << AmpSq0.u << endl ;
	
	double integrand ;
	if (H.IsParton_u() == true) integrand = H.AmpSq().u ;	
	else if (H.IsParton_d() == true) integrand = H.AmpSq().d ;	
	
	return integrand*Jac ;
}

partons<double> xsection::Hard()
{
	double mu2 ;
	if (Nloop > 0) mu2 = mu_h2 ;
	else mu2 = mu_r2 ;

	partons<bool> flag ;	
	double result, error ;
	hard H ;
	if ((Nloop == 1 && part == 1) || Nloop == 2)
	H.set_para (mu2,  M2, 0) ;
	else H.set_para (mu2,  M2, Nloop) ;
	partons<double> AmpSq ;
	
	//printf("M = % .3f , mu = % .3f \n" , sqrt(M2), sqrt(mu2)) ;
	
	H.select_u() ;
	NIntegrate(0., 1., &H, integrand_cosT, result, error, 0);
	AmpSq.u = result ;
	
	H.select_d() ;
	NIntegrate(0., 1., &H, integrand_cosT, result, error, 0);
	AmpSq.d = result ;
	
	return AmpSq;
}

/*
//==========================================================================
partons<double> xsection::AmpSq()
{
	double mu2 ;
	if (Nloop == 2) mu2 = mu_RG2 ;
	else mu2 = mu_r2 ;
	
	partons<double> ampLO = AmpLO() ;
	if (Nloop == 0) return ampLO ;
	else if (Nloop == 1) {
		partons<double> AmpNLO ;
		AmpNLO.u = ampLO.u*real(2.*WilsonC(mu2).u-1.) ; 
		AmpNLO.d = ampLO.d*real(2.*WilsonC(mu2).d-1.) ; 
		AmpNLO.g = ampLO.g*real(2.*WilsonC(mu2).g-1.) ;
		return AmpNLO ;
	} 
	else if (Nloop == 2){
		if (Matching == false) {
			mu_h2 = -M2 ;
			double anom = Col_anomaly() ;
			//printf("anomaly = %6.3f \n",anom) ;
			dcomplex run = WilsonC_evolve(mu_h2, mu2, M2) ;
			//printf("Wilson running = %6.3f \n", abs(run)) ;
			//abort() ; 
			partons<double> AmpNLL ;
			AmpNLL.u = ampLO.u*norm(run*WilsonC(mu_h2).u)*anom ; 
			AmpNLL.d = ampLO.d*norm(run*WilsonC(mu_h2).d)*anom ; 
			AmpNLL.g = ampLO.g*norm(run*WilsonC(mu_h2).g)*anom ;
			return AmpNLL ;
			}
		else {
			double anom = Col_anomaly() ; 
			partons<double> AmpNLL ;
			AmpNLL.u = ampLO.u*real(2.*WilsonC(mu2).u +anom -2.) ; 
			AmpNLL.d = ampLO.d*real(2.*WilsonC(mu2).d +anom -2.) ; 
			AmpNLL.g = ampLO.g*real(2.*WilsonC(mu2).g -1.) ;
			return AmpNLL ;
		}
	}
}

//==========================================================================
double xsection::Col_anomaly()
{
	if (Matching == false) {
		double fqq = Fqq(pT_veto, sqrt(mu_RG2), R, 2) ;
		return pow( M/pT_veto, -2*fqq) ;
	}
	else {
		double fqq = Fqq(pT_veto, sqrt(mu_RG2), R, 1) ;
		return (1. - 2*fqq*log( M/pT_veto)) ;
	}
}

//==========================================================================
//  Function calculates pre-factor in the xsection formula
double xsection::factor()
{
	double factor = 1.0/NC/NC * pbarn
				 * betaW/16/pi/M/ECM/ECM 
				 * 2*GF*GF*MW2*MW2 ;
	return factor ;
}

//==========================================================================
//  Integrand of the differential cross-section in M, Y and cos(theta)
// The limits for x are taken from 0 to 1 and the appropriate Jacobian included
double integrand_xs(double* x, size_t Ndim, void* para_ptr)
{ 
	xsection xs = *(xsection*) para_ptr ;
	double Jac = 1.0 ;  // Take care of Jacobian in the integrals
	
	//Integration over the invariant mass of the WW pair
	double ECM = xs.get_ECM(), ECM2 = ECM*ECM ;
	double M_min = 2*MW, M_max = ECM ;
	double M = M_min + x[0]*(M_max - M_min) ;
	double M2 = M*M ;
	Jac = Jac*(M_max - M_min) ;
	
	//Integration over the rapidity of the WW pair
	double Y_min = -log(ECM/M) ;
	double Y_max = log(ECM/M) ;
	double Y = Y_min + x[1]*(Y_max - Y_min) ;
	Jac = Jac*(Y_max - Y_min) ;
	
	// Integration over cos(theta) 
	// theta is the angle between p3 and the x-axis in CoM frame
	double cosT_min = -1.0 ;
	double cosT_max = 1.0 ;
	double cosT = cosT_min + x[2]*(cosT_max - cosT_min) ;
	Jac = Jac*(cosT_max - cosT_min) ; 
	
	
	//set the phase space for the xsection :
	xs.set_PhaseSpace(M, Y, cosT) ;

	
	partons<double> lum = xs.parton_lum() ;
	partons<double> AmpSq = xs.AmpSq() ;
	double factor = xs.factor() * Jac ;
	
	double integrand ;
	//if (xs.get_Nloop() == 2 && xs.isMatching() == true)
	//integrand = factor * 
	//	 		(lum.u * xs.AmpLO().u + lum.d * xs.AmpLO().d)  ;
	//else
	integrand = factor * 
		 		(lum.u * AmpSq.u + lum.d * AmpSq.d)  ;
	
	//cout << WilsonC.u << " " << WilsonC.d << " " << WilsonC.g << endl ;
	//cout << factor << " " << lum.u << " " << WilsonC.u << endl ;
	//cout << lum.u << endl ;
	//factor * lum.u * WilsonC.u ; 
	//cout << " M , Y, cosT : " << M << " " << Y << " " << cosT << endl ;
	
//	double integrand = 8.0*x[0]*x[1]*x[2];
	return integrand ;
}
*/

//==========================================================================
//  Function calculates pre-factor in the xsection formula
double xsection::factor()
{

	double factor = 1.0/NC/NC*pbarn/16/pi/ECM/ECM 
				 * 2*GF*GF*MW2*MW2 ;
	return factor ;
}


//==========================================================================
double xsection::Col_anomaly()
{
	if (Nloop >= 2) {
		double fqq, Nu_fac ;
		Nu_fac = 1. ;
		double mu = sqrt(mu_RG2) ;
		if (Nloop == 4) 
		fqq = Fqq(pT_veto, mu, R, 1) ;
		else 
		fqq = Fqq(pT_veto, mu, R, Nloop-1) ;
		double hf;
		if (Nloop == 2 || Nloop == 4)  hf = 0.;
		else  hf = hF(pT_veto, mu) ; 
		//cout << "fqq : " << fqq << "     Rap : " <<  pow(NuScaleFactor, -2*fqq) << endl ;
		double as = alpha_s(mu*mu)/4/pi ;
		double LT = 2.*log(mu/pT_veto) ;
		if (Nloop == 3 ) 
			Nu_fac = 1. + as*G0*LT*log(NuScaleFactor) ;
		return pow(sqrt(NuScaleFactor) * M/mu, -2.*fqq ) * exp(2*hf) * Nu_fac ;
		//return pow(M/mu, -2*fqq) * exp(2*hf) ;
	}
	else if (Nloop == 1) {
		double mu = sqrt(mu_RG2) ;
		double fqq = Fqq(pT_veto, mu, R, 1) ;
		double hf = hF(pT_veto, mu) ; 
		//return (- 2*fqq*log( M/pT_veto) + 2*hf) ;
		return (- 2*fqq*log( M/mu) + 2*hf) ;
	}
	else return 1. ;
}

double integrand_diffM(double x, void* para_ptr)
{
	xsection xs = *(xsection*) para_ptr ;
	double Jac = 1.0 ;  // Take care of Jacobian in the integrals
	
	//Integration over the invariant mass of the WW pair
	double ECM = xs.get_ECM(), ECM2 = ECM*ECM ;
	double M_min = 2*MW, M_max = ECM ;
	double M = M_min + x*(M_max - M_min) ;
	double M2 = M*M ;
	Jac = Jac*(M_max - M_min) ;
	
	double betaW = sqrt(1.0 -4*MW2/M2) ;
	//printf("M = % .3f  \n" , M) ;
	xs.M = M ; xs.M2 = M*M ;
	xs.set_scales() ;
	if (xs.dynamic_scale) 
		xs.set_muh2(M2*xs.scaleFactor) ;
	
	double integrand ;
	
	if (xs.Nloop == 1) {
		xs.part = 0 ;
		partons<double> Lum_0, Hard_0, Lum_1, Hard_1 ; 
		Lum_0 = xs.Lum() ; Hard_1 = xs.Hard() ;
		xs.part = 1 ;
		Lum_1 = xs.Lum() ; Hard_0 = xs.Hard() ;
		double col_anom = xs.Col_anomaly() ;
//		double run = 2.*real(WilsonC_evolve(xs.mu_h2, xs.mu_RG2, 
//							xs.M2, 0)) ;
		double muff2 = xs.mu_RG2 ;
		double run = -CF*alpha_s(muff2)/2/pi*(pow(log(M2/muff2),2.)
							- 3.*log(M2/muff2)) ;
		double runPi = -CF*alpha_s(muff2)/2/pi*(pow(log(M2/muff2),2.)
							-pi*pi - 3.*log(M2/muff2)) ;
		double hf = 2.*hF(xs.pT_veto, sqrt(muff2)) ;
		integrand = (Lum_0.u*Hard_1.u + Lum_0.d*Hard_1.d) 
				 + (Lum_1.u*Hard_0.u + Lum_1.d*Hard_0.d) 
		 		+ (Lum_0.u*Hard_0.u + Lum_0.d*Hard_0.d)
					*(col_anom) ;

	if (xs.is_diff) {
	cout << Hard_0.u << "    " << Hard_0.d << "    "  << Hard_1.u 
	<< "    " << Hard_1.d << "    " << Lum_0.u <<  "    " << Lum_0.d 
	<< "    " << Lum_1.u <<  "    " << Lum_1.d << "    " 
	<< col_anom <<  "    " << run << "    " << runPi << "    " 
	<< hf << "    " ;
	}
	}
	else if (xs.Nloop >= 2) {
		partons<double> Lum, Hard ;
		double col_anom = xs.Col_anomaly() ;
		double run = norm(WilsonC_evolve(xs.mu_h2, xs.mu_RG2, 
							xs.M2, xs.Nloop-1)) ;
		Lum = xs.Lum() ; Hard = xs.Hard() ;
		integrand = (Lum.u*Hard.u + Lum.d*Hard.d)
					* col_anom * run ;
	if (xs.is_diff) {
	cout << Hard.u*run << "    " << Hard.d*run << "    " << Lum.u << 
	"    " << Lum.d << "    " << col_anom <<  "    " ;
	}
	}
	else if (xs.Nloop == 0) {
		partons<double> Lum, Hard ;
		Lum = xs.Lum() ; Hard = xs.Hard() ;
		integrand = (Lum.u*Hard.u + Lum.d*Hard.d) ;
	if (xs.is_diff) {
	cout << Hard.u << "    " << Hard.d << "    " << Lum.u << 
	"    " << Lum.d << "    "  ;
	}
	}
	//cout << Lum.u << endl ;
	//abort() ;
	//partons<double> Hard = xs.Hard() ;
	//printf("M = % .2f , Lum = % .5f , Hard = % .3f \n" , M, Lum.u , Hard.u) ;
	if (xs.is_diff) {
	double factor = 1.0/NC/NC*pbarn*1.0/16/pi/ECM/ECM *2*GF*GF*MW2*MW2*
	1000.*betaW/M ;
	cout << factor << "    "  ;
	}	
	return Jac*betaW/M*integrand ;
}

double evaluate_xs_diffM(double M, xsection* para_ptr) 
{
	xsection xs = *para_ptr ;
	double ECM = xs.get_ECM() ; 
	double factor = 1.0/NC/NC*pbarn
					*1.0/16/pi/ECM/ECM 
				    *2*GF*GF*MW2*MW2 ;
	double M_min = 2*MW, M_max = ECM ;
	double x = (M - M_min)/(M_max - M_min) ;
	return factor*1000.*integrand_diffM(x, &xs)/(M_max-M_min);
}

// Function to put together the final cross-section
double evaluate_xs(xsection* para_ptr) 
{
	
	double result, error ;
	xsection xs = *para_ptr ;
	double ECM = xs.get_ECM() ; 
	double factor = 1.0/NC/NC*pbarn
					*1.0/16/pi/ECM/ECM 
				    *2*GF*GF*MW2*MW2 ;
	NIntegrate(0., 1., &xs, integrand_diffM, result, error, 0);
	//printf("factor = % .5f , xsec = % .3f \n" , factor , result) ;
	//	cout << "factor = " << factor << endl ;
	return result*factor ;
}

//==========================================================================
int main (void)
{

	timer t ;

/*	
	xsection xs(2,2) ;
	xs.set_R(0.5) ;
	xs.set_ECM(8000.) ;
	xs.set_dynamic(true) ;
	xs.set_muf2(MW2) ; xs.set_mur2(MW2) ;
	xs.set_muh2(MW2) ; 
	
	ofstream myfile;
  	double pTveto, S0, S1, S2, S3, S4 ;
	for (int i = 1 ; i < 17 ; i++) {
		t.start() ;
		pTveto = 4.+2.*double(i) ;
		myfile.open ("Results/NNLL.dat", ios::out | ios::app);
		cout << "===========================\n" ;
		cout << "Processing for pT_veto = " << pTveto << endl;
		myfile.width(8) ; myfile.precision(3) ; 
		myfile << pTveto ; myfile.precision(4) ;
		xs.set_scaleFactor(-1.) ;
		xs.set_pTveto(pTveto); 
		myfile.width(8) ;
		xs.set_muRG2(pTveto*pTveto); S0 = evaluate_xs(&xs);
		printf("central : %8.3f \n" , S0) ; myfile << S0 ;
		myfile.width(8) ;
		xs.set_muRG2(pTveto*pTveto/4.); S1 = evaluate_xs(&xs);
		printf("pTveto down : %8.3f \n" , S1) ; myfile << S1;
		myfile.width(8) ;
		xs.set_muRG2(pTveto*pTveto*4.); S2 = evaluate_xs(&xs);
		printf("pTveto up : %8.3f \n" , S2) ; myfile << S2 ;
		myfile.width(8) ;
		xs.set_muRG2(pTveto*pTveto);
		xs.set_scaleFactor(-1/4.); S3 = evaluate_xs(&xs);
		printf("hard down : %8.3f \n" , S3) ; myfile << S3;
		myfile.width(8) ;
		xs.set_scaleFactor(-4.); S4 = evaluate_xs(&xs);
		printf("hard up : %8.3f \n" , S4) ; myfile << S4 ;
		myfile << "\n" ;
		t.stop() ;
		cout << "Time taken : " << t.duration() << " secs" << endl;
		myfile.close();
	}
 */



	vector<string> input_list ;
	input_list.clear() ;
	string input ;
	int comment_pos ;
	bool eof = false ;
	while(!eof) {
		getline(cin, input) ;
		comment_pos = input.find_first_of("#") ;
		if (comment_pos == 0) continue ;
		else if (comment_pos > 0) input.resize(comment_pos) ;
		input_list.push_back(input) ;
		if (input == "" || input[0] == ' ') eof = true ;
	}	
	
	int i = 0 ;
	int is_diff = atof(input_list[i++].c_str()) ;
	double ECM = atof(input_list[i++].c_str()) ; 
	int N = atof(input_list[i++].c_str()) ;
	int N_pdf = atof(input_list[i++].c_str()) ;
	int N_f = atof(input_list[i++].c_str()) ;
	double R = atof(input_list[i++].c_str()) ;
	double mu_r = atof(input_list[i++].c_str()) ;
	double mu_f = atof(input_list[i++].c_str()) ;
	double pTveto = atof(input_list[i++].c_str()) ;
	double mu_RG = atof(input_list[i++].c_str()) ;
	int dynamic = atof(input_list[i++].c_str()) ;
	double mu_h = atof(input_list[i++].c_str()) ;
	double ScaleFac = atof(input_list[i++].c_str()) ;
	double NuScaleFac = atof(input_list[i++].c_str()) ;

			
	xsection xs(N,N_pdf) ;
	xs.set_ECM(ECM) ;
	xs.set_R(R) ;
	xs.set_pTveto(pTveto); 
	xs.set_dynamic(dynamic) ;
	//xs.set_scaleFactor(ScaleFac) ;
	xs.set_scaleFactor(ScaleFac, NuScaleFac) ;
	xs.set_muf2(mu_f*mu_f) ; xs.set_mur2(mu_r*mu_r) ;
	xs.set_muh2(ScaleFac*mu_h*mu_h) ; 
	xs.set_muRG2(mu_RG*mu_RG) ;
	xs.set_diff(is_diff) ;
	if (N_f == 5) xs.incl_b(true) ;
	else xs.incl_b(false) ;
	


	if (is_diff) {
//	cout << "160.4   0.0" << endl ; 			    
	for (double M = 161.; M <= 175.; M = M+1.){
		cout << M << "  " ;
		cout << evaluate_xs_diffM(M, &xs) << endl;
	}
	for (double M = 180.; M <= 230.; M = M+5.){
		cout << M << "  " ;
		cout << evaluate_xs_diffM(M, &xs) << endl;
	}
	for (double M = 240.; M <= 610.; M = M+10.){
		cout << M << "  " ;
		cout << evaluate_xs_diffM(M, &xs) << endl;
	}
	}
	else cout << evaluate_xs(&xs) << endl ;



		
	
/*
	xsection xs(1,2) ;
	xs.set_R(0.5) ;
	xs.set_dynamic(false) ;
	xs.set_muf2(MW2) ; xs.set_mur2(MW2) ;
	xs.set_muh2(MW2) ;

	ofstream myfile;
	double pTveto, S0, S1, S2 ;
	for (int i = 3 ; i < 17 ; i++) {
		t.start() ;
		pTveto = 4.+2.*double(i) ;
		xs.set_pTveto(pTveto);
		myfile.open ("Results/NNLL_expansion.dat", ios::out | ios::app);
		cout << "===========================\n" ;
		cout << "Processing for pT_veto = " << pTveto << endl;
		myfile.width(8) ; myfile.precision(3) ; 
		myfile << pTveto ; myfile.precision(4) ;
		xs.set_muh2(MW2) ; xs.set_muRG2(MW2);  
		myfile.width(8) ;
		S0 = evaluate_xs(&xs);
		printf("central : %8.3f \n" , S0) ; myfile << S0 ;
		myfile.width(8) ;
		xs.set_muh2(MW2/4.) ; xs.set_muRG2(MW2/4.); 
		S1 = evaluate_xs(&xs);
		printf("hard down : %8.3f \n" , S1) ; myfile << S1;
		myfile.width(8) ;
		xs.set_muh2(MW2*4.) ; xs.set_muRG2(MW2*4.); 
		S2 = evaluate_xs(&xs);
		printf("hard up : %8.3f \n" , S2) ; myfile << S2 ;
		myfile << "\n" ;
		t.stop() ;
		cout << "Time taken : " << t.duration() << " secs" << endl;
		myfile.close();
	}
*/

/*	
	double result, error; 
	params_lum para_lum ;
	para_lum.mu = MW ;
	double tau = 0.1 ;
	para_lum.tau = tau ;
		para_lum.partons_ID = 0 ;
		NIntegrate(0., 1., &para_lum, 
				integrand_lum, result, error, 0);
	cout << result << endl;
*/
	//cout <<  evaluate_xs(&xs) << endl;
/*
	for (int i = 0; i < 11; i++) {
		double pTveto = pow(10.,i/5.+1) ;
		xs.set_pTveto(pTveto) ;
		cout << pTveto << "  " << evaluate_xs(&xs) << endl;
	}
*/
/*
	for (int i = 2; i < 13; i++) {
		double pTveto = 2.5*i ;
		xs.set_pTveto(pTveto) ;
		cout << pTveto << "  " << evaluate_xs(&xs) << endl;
	}
	for (int i = 6; i < 15; i++) {
		double pTveto = 5. + 5.*i ;
		xs.set_pTveto(pTveto) ;
		cout << pTveto << "  " << evaluate_xs(&xs) << endl;
	}
*/
	//cout << hF(MW, MW) << endl ;
	//cout << norm(WilsonC_evolve(MW2, MW2, 10*MW2)) << endl ;
	//cout << Fqq(MW, MW, 0.4, 2) << endl ;
	//xs.set_Match(true) ;
/*
	double result, error ;
	result = 0.0 ; error = 0.0 ;
//	double a[3] = {1.0 , 2.0, 3.0} ;
	double lower_limit[3] = {0.0, 0.0, 0.0};
	double upper_limit[3] = {1.0, 1.0, 1.0};
	void* params_ptr = &xs ;
	t.start() ;
	NIntegrate(lower_limit,upper_limit, 3, params_ptr, 
			integrand_xs, result, error, 3000, 1);
	//cout << " Result = " << result << endl ;
	//cout << " MW = " << MW << endl ;
	//xs.set_PhaseSpace(200., 0.5, 0.2) ;
	//xs.set_scales() ;
	//cout << "full anom " << xs.Col_anomaly() << endl ;
	//xs.set_Match(true) ;
	//xs.set_scales() ;
	//cout << "approx anom " << xs.Col_anomaly() << endl ;
	double S = 80000. , T = -4000. ;
	double U = 2* MW2 - S - T ;
	//xs.set_Mandelstam(S,T,U) ;
	//cout << xs.WilsonC(MW2).u << endl ;
	//cout << " wilson : " << endl ;
	//cout << F0_SS(S,T,U) << endl ;
	//cout << F0_ST(S,T,U) << endl ;
	//cout << F0_TT(S,T,U) << endl ;
	//cout << F1_SS(S,T,U,MW) << endl ;
	//cout << F1_ST(S,T,U,MW) << endl ;
	//cout << F1_TT(S,T,U,MW) << endl ;
	//cout << " couplings : " << endl ;
	//cout << C_SS_u(S,T,U) << endl ;
	//cout << C_ST_u(S,T,U) << endl ;
	//cout << C_SS_d(S,T,U) << endl ;
	//cout << C_ST_d(S,T,U) << endl ;
	//cout << "final : " << endl ;
	//double Mlo = 1./16/SW2/SW2/NC/NC *F0_ST(S,T,U)*C_ST_u(S,T,U) ;
	//double Mlo = 1./16/SW2/SW2/NC/NC *F0_TT(S,T,U) ;
	//double Mlo = 1./16/SW2/SW2/NC/NC *(F0_SS(S,T,U)*C_SS_u(S,T,U)
	 //+ 2*F0_ST(S,T,U)*C_ST_u(S,T,U) + F0_TT(S,T,U));
	//cout << "LO " << Mlo <<endl;
	//dcomplex Cnlo = 4.*F1_TT(S,T,U,MW2,0)/F0_TT(S,T,U) ;
	//dcomplex Cnlo = 4.*(F1_SS(S,T,U,MW2,-2)*C_SS_u(S)
	 //+ F1_ST(S,T,U,MW2,-2)*C_ST_u(S) + F1_TT(S,T,U,MW2,-2))/
	 //(F0_SS(S,T,U)*C_SS_u(S)
	 //+ 2*F0_ST(S,T,U)*C_ST_u(S) + F0_TT(S,T,U)) ;
	//cout << "NLO " << Cnlo <<endl;
	//printf("ratio alphas % .6f + i % .6f \n", real(alphas(-MW2)/alphas(MW2)), 
	//imag(alphas(-MW2)/alphas(MW2))) ;
	//printf("alphas % .3f , % .3f \n" , real(alphas(MW2)) , alpha_s(MZ2));
	//cout << "dilog " << gsl_sf_dilog(-1000.) << endl ;
	//cout << 1./16/SW2/SW2/NC/NC *(F0_SS(S,T,U)*C_SS_u(S,T,U) 
	//		+ F0_TT(S,T,U) + 2*F0_ST(S,T,U)*C_ST_u(S,T,U)) <<endl;
	//xs.set_PhaseSpace(200., 0.5, 0.2) ;
	//xs.set_scales() ;
	//cout << "alphas " << alpha_s(MW) <<endl;
	//cout << "Phase space " << xs.factor() << endl ;
	//cout << "Parton lum " << xs.parton_lum().d << endl ;
*/

	return 0; 
	
}
