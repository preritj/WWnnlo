#ifndef WW_H
#define WW_H

#include <const.h>
#include <process.h>

using namespace smPara;

class pp2WW : public process{
	private :
		flav T3, Q, PosUp ;
		flav CQtt ;	
		double coup ;
	public :
		pp2WW() ;
		flav CQst(), CQss() ;
		double F0(), J0(), K0() ;
		double F1(), J1(), K1() ;
		void calc_lo(), calc_nlo() ; //Calculate H0 and H1 
};
#endif
