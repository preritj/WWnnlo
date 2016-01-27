#ifndef XSEC_H
#define XSEC_H

#include <vector>
#include <process.h>
#include <Hard.h>
#include <Beam.h>
#include <const.h>

using namespace beamPara;

class Xsection{
	private : 
		double mu, as;
		process& proc;
		Hard hard;
		Beam beam;
	public :
		Xsection(process& proc_) : 
			proc(proc_), hard(proc_) {
				set_mu(pTveto) ;	
			}
		double factor(double);
		void set_mu(double) ;
		flav lo(double), lo(double, double) ;
		vector<flav> nlo(double), nlo(double, double) ;
};
#endif
