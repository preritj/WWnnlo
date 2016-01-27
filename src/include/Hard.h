#ifndef HARD_H
#define HARD_H

#include "vector"
#include "flavor.h"
#include "process.h"

struct hPara{
	int order;			// loop-order LO:0, NLO:1
	char flavor;		// 'u' for up-type, 'd' for down-type
	double M ; 			// invariant mass of final state excluding radiation
	process& proc;		// reference to process class object
	hPara(process& proc_) : proc(proc_) {}
};


class Hard{
	private : 
		double  mu ;
		process& proc;
	public : 
		Hard(process& proc_) : proc(proc_) {}
		void set_mu(double);
		flav H0(), H1reg();
		double  H1mu(double, double), H1mu(double) ;
		flav lo(double);
		vector<flav>  nlo(double, double), nlo(double) ;
};

double integrand_cosT(double , void* ) ;
#endif
