#ifndef WW_H
#define WW_H

#include <input.h>
#include <Hard.h>

class WW : public Hard{
	private :
		double s, t, u ;
	public :
		WW(InputPara input) : Hard(input) {}
		double F0(); double J0(); double K0() ;
		void set_Mandelstam(double, double);
		void set_Mandelstam(double, double, double);
};
#endif
