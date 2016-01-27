// cos(theta) Integration
//-------------------------------------------------------------
#include <iostream> 
#include <cstdlib> 
#include <math.h>
#include <string>
#include <vector>
#include <algorithm>
#include <const.h>
#include <misc.h>
#include <Hard.h>
#include <process.h>

using namespace std;
using namespace smPara;

double integrand_cosT(double x, void* para_ptr){
	double integrand = 1. ;
	double Jac = 1.0 ;  // Take care of Jacobian in the integrals
	hPara para = *(hPara*) para_ptr ;
	
	double cosT_min = -1.;
	double cosT_max =  1. ;
	double cosT = cosT_min + x*(cosT_max - cosT_min) ;
	Jac = Jac * (cosT_max - cosT_min) ;

	int n = para.order ;
	char flavor = para.flavor ;
	process& proc = para.proc ;
	double M1 = proc.ext_masses()[0]; 
	double M2 = proc.ext_masses()[1]; 
	double s = pow(para.M, 2);
	double dr = (M1*M1 - M2*M2)/s ;
	double beta = sqrt(1.0 + dr*dr - 2/s*(M1*M1+M2*M2) ) ;
	double t  = (M1*M1 + M2*M2)/2 - s/2* (1 - beta * cosT);
	proc.set_Mandelstam(s,t) ;
	switch(n){
	case 0 :
		if (flavor == 'u') integrand = proc.Hlo().u ;
		else if (flavor == 'd') integrand = proc.Hlo().d ;
		else {cout << "Invalid flavor in hard parameters. Aborting..." << endl;
			abort() ;}
		break;
	case 1 :
		if (flavor == 'u') integrand = proc.Hnlo_reg().u ;
		else if (flavor == 'd') integrand = proc.Hnlo_reg().d ;
		else {cout << "Invalid flavor in hard parameters. Aborting..." << endl;
			abort() ;}
		break;
	default :
		cout << "Invalid loop order in hard parameters. Aborting..." << endl;
		abort() ;
	}
	return integrand*Jac ;
}
