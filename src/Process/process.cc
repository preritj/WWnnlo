#include <const.h>
#include <math.h>
#include <process.h>
#include <iostream>

using namespace std;
using namespace smPara;

// set Mandelstam variables
void process::set_Mandelstam(double s_, double t_){
	double u_ = M1*M1 + M2*M2 -s_ - t_ ;
	set_Mandelstam(s_,t_,u_) ;
}

void process::set_Mandelstam(double s_, double t_, double u_){
	s=s_; t=t_; u=u_; 
	double dr = (M1*M1 - M2*M2)/s ;
	beta=sqrt(1.+dr*dr-2/s*(M1*M1+M2*M2));
}

//return LO hard functions
flav process::Hlo(double s_, double t_, double u_){
	set_Mandelstam(s_,t_, u_);
	return Hlo() ;
}

flav process::Hlo(double s_, double t_){
	set_Mandelstam(s_,t_);
	return Hlo() ;
}

//return NLO hard functions
flav process::Hnlo_reg(double s_, double t_, double u_){
	set_Mandelstam(s_,t_, u_);
	return Hnlo_reg() ;
}

flav process::Hnlo_reg(double s_, double t_){
	set_Mandelstam(s_,t_);
	return Hnlo_reg() ;
}

vector<double> process::ext_masses(){
	vector<double> m ; m.clear() ;
	m.push_back(M1) ;
	m.push_back(M2) ;
	return m;
}
