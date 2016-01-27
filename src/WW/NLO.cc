#include <const.h>
#include <WW.h>
#include <math.h>
#include <gsl/gsl_sf_dilog.h>

using namespace constants;
using namespace smPara;

double pp2WW::F1(){
	return 
	(64*t*(s + t))/pow(MW2,2) - (128*(2*s + t))/MW2 - (4*pow(s + 
	4*t,2))/(pow(beta,2)*s*t) + (4*(72*pow(MW2,2) - 140*MW2*t + 73*s*t + 
	80*pow(t,2)))/pow(t,2) + (32*pow(pi,2)*(-((t*(s + t))/pow(MW2,2)) + 
	(2*(2*s + t))/MW2 - (-4*MW2 + 2*s + 3*t)/t))/3. + ((6*pow(s + 
	4*t,2))/(pow(beta,4)*s*t) + (8*(12*pow(MW2,2) - 19*MW2*t + 8*s*t + 
	6*pow(t,2)))/pow(t,2) - (-26*pow(s,2) - 128*s*t + 
	32*pow(t,2))/(pow(beta,2)*s*t))*log(s/MW2) - ((128*s)/(-MW2 + t) + 
	(32*(-3*pow(MW2,2) - 3*s*t + pow(t,2)))/pow(t,2))*log(-(t/MW2)) - 
	(64*(-MW2 + t)*((2*pow(MW2,2))/pow(t,2) - u/t)*(-pow(pi,2)/2. + 
	pow(log(MW2/s),2)/2. - pow(log(-(t/s)),2)/2.))/(MW2 - t) + 
	(((-3*pow(s + 4*t,2))/(2.*pow(beta,4)*t) + (2*(-3*pow(s,2) - 14*s*t + 
	8*pow(t,2)))/(pow(beta,2)*t) + (-48*pow(MW2,2) + 72*MW2*s - 
	49*pow(s,2) + 16*t*(4*MW2 - u))/(2.*t))*(pow(pi,2)/3. + pow(log((1 - 
	beta)/(1 + beta)),2) + 4*gsl_sf_dilog((-1 + beta)/(1 + 
	beta))))/(beta*s) + (32*((2*pow(MW2,2))/t - u)*(-4*log((MW2 - 
	t)/MW2)*log(-(t/MW2)) + 2*pow(log(-(t/MW2)),2) - 
	4*gsl_sf_dilog(t/MW2)))/t;
}

double pp2WW::J1(){
	return
	((32*pow(s,2))/(-MW2 + t) + (64*s*t*(s + t))/pow(MW2,2) - 
	(16*(17*pow(MW2,2) + 34*MW2*s - 26*MW2*t - 21*s*t + pow(t,2)))/t - 
	(128*(2*pow(s,2) + 2*s*t + pow(t,2)))/MW2 + (32*pow(pi,2)*(-4*s - t - 
	(s*t*(s + t))/pow(MW2,2) - (2*MW2*(-MW2 - 2*s + 2*t))/t + 
	(2*(2*pow(s,2) + 2*s*t + pow(t,2)))/MW2))/3. + ((48*MW2*(-MW2 - 2*s + 
	2*t))/t - 16*(-2*s + 3*t) + (16*(s + 4*t))/pow(beta,2))*log(s/MW2) + 
	((-48*MW2*(MW2 + 2*s))/t - (32*pow(s,2)*t)/pow(-MW2 + t,2) + 
	16*(2*MW2 - 5*s + t) + (64*s*(s + 2*t))/(-MW2 + t))*log(-(t/MW2)) + 
	(32*(-MW2 + t)*(-2*s + (2*MW2*(MW2 + 2*s))/t - u)*(-pow(pi,2)/2. + 
	pow(log(MW2/s),2)/2. - pow(log(-(t/s)),2)/2.))/(MW2 - t) + 
	((32*pow(MW2,2) - 12*pow(s,2) + 32*s*t - 16*MW2*(7*s + 2*t) - (4*s*(s 
	+ 4*t))/pow(beta,2))*(pow(pi,2)/3. + pow(log((1 - beta)/(1 + 
	beta)),2) + 4*gsl_sf_dilog((-1 + beta)/(1 + beta))))/(beta*s) + 
	(16*(-2*MW2*(MW2 + 2*s) + t*(2*s + u))*(-4*log((MW2 - 
	t)/MW2)*log(-(t/MW2)) + 2*pow(log(-(t/MW2)),2) - 
	4*gsl_sf_dilog(t/MW2)))/t)/s;
}

double pp2WW::K1(){
	return 
	16*(2 - pow(pi,2)/3.)*(12*pow(MW2,2) - 4*MW2*s + 17*pow(s,2) - 
	24*MW2*t + 20*s*t + 12*pow(t,2) + (pow(s,2)*t*(s + t))/pow(MW2,2) - 
	(2*s*(2*pow(s,2) + 3*s*t + 2*pow(t,2)))/MW2)/s/s;
}
