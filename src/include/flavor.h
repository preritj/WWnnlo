#ifndef FLAV_H
#define FLAV_H

#include <math.h>

struct flav{
	double u, d;	 // up and down type
};

inline flav operator+(const flav& f1, const flav& f2){
	flav f;
	f.u = f1.u + f2.u;
	f.d = f1.d + f2.d;
	return f;
}

inline flav operator-(const flav& f1, const flav& f2){
	flav f;
	f.u = f1.u - f2.u;
	f.d = f1.d - f2.d;
	return f;
}

inline flav operator*(const flav& f1, const flav& f2){
	flav f;
	f.u = f1.u * f2.u;
	f.d = f1.d * f2.d;
	return f;
}

inline flav operator*(const flav& f1, double factor){
	flav f;
	f.u = f1.u*factor;
	f.d = f1.d*factor;
	return f;
}

inline flav operator*(double factor, const flav& f1){
	flav f;
	f.u = f1.u*factor;
	f.d = f1.d*factor;
	return f;
}

inline flav pow(const flav& f1, double power){
	flav f;
	f.u = pow(f1.u, power);
	f.d = pow(f1.d, power);
	return f;
}
#endif
