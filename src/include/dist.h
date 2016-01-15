#ifndef DIST_H
#define DIST_H

struct Dist{
	friend Dist operator+(const Dist&, const Dist&) ;
	friend Dist operator*(const Dist&, double) ;
	Dist(){ None=0.; Plus=0.; Delta=0.;}
	double None, Plus, Delta;
};

Dist Iqfromq(double , double );
Dist Iqfromg(double , double );
double PlusRem(double);

#endif
