#include <math.h>
#include <PDF.h>

// Returns PDF in f
void evolvePDF(double x, double Q, double* f)
{
	evolvepdf_(x , Q , f) ;
}

//running QCD coupling
double alpha_s (double mu)
{
	return alphaspdf_(mu) ;
}
