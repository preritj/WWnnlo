#include <WW.h>
#include <iostream>

void WW::set_Mandelstam(double s_, double t_){
	double u_ = 2*MW2 -s_ - t_ ;
	set_Mandelstam(s_,t_,u_) ;
}

void WW::set_Mandelstam(double s_, double t_, double u_){
	s=s_; t=t_; u=u_;
}
