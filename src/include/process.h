#ifndef PROCESS_H
#define PROCESS_H

#include <flavor.h>
#include <vector>

using namespace std;

class process{
	protected :
		double s, t, u, M1, M2, beta;
		flav H0, H1;
	public :
		void set_Mandelstam(double, double);
		void set_Mandelstam(double, double, double);
		vector<double> ext_masses(); 
		flav Hlo() {calc_lo(); return H0;};
		flav Hnlo_reg() {calc_nlo(); return H1;};
		flav Hlo(double, double), Hnlo_reg(double, double);
		flav Hlo(double, double, double), Hnlo_reg(double, double, double);
		virtual void calc_lo(){H0.u = 0.; H0.d = 0.;};
		virtual void calc_nlo(){H1.u = 0.; H1.d = 0.;};
};
#endif
