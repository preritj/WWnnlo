#include <iostream> 
#include <fstream> 
#include <complex>
#include <stdlib.h> 
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <vector>

#include <misc.h>
#include <PDF.h>
#include <Beam.h>
#include <input.h>
#include <WW.h>
#include <Hard.h>
#include <Xsection.h>
#include <const.h>

using namespace std;
using namespace smPara;
using namespace beamPara;


int main (void)
{

	timer t ;
	t.start();


	// User code goes here
	//------------------------------------------------------- 
	double M = 500. ;
	Beam B ;
	cout << "Invariant mass of W pair, M =  " << M << endl ;
	cout << "Beam-up LO    = " << B.lo(500).u << endl;  
	cout << "Beam-down LO  = " << B.lo(500).d << endl;  
	cout << "Beam-up NLO   = " << B.nlo(500).u << endl;  
	cout << "Beam-down NLO = " << B.nlo(500).d << endl;  
	pp2WW ww;
	ww.set_Mandelstam(24000.,-1200.);
	Hard hard(ww);
	vector<flav> H = hard.nlo(M, 25.); 
	cout << "H0(M) = " << H[0].u << endl;
	cout << "H1(M) = " << H[1].u << endl;
	Xsection xs(ww);
	xs.set_mu(pTveto);
	flav nlo = xs.nlo(M)[1];
	cout << " NLO(u) diff-xsection (pb) = " << nlo.u << endl;
	cout << " NLO(d) diff-xsection (pb) = " << nlo.d << endl;
	// User code ends here
	//------------------------------------------------------- 
	t.stop();
	cout << "Time used : " << t.duration() << " secs" << endl;
	return 0;
}	
