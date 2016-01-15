#include <iostream> 
#include <fstream> 
#include <complex>
#include <stdlib.h> 
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <vector>

#include <misc.h>
#include <functions.h>
#include <Beam.h>
#include <input.h>

using namespace std;


int main (void)
{

	timer t ;
	t.start();


	// User code goes here
	//------------------------------------------------------- 
	double M = 500. ;
	ReadInput input;
	Beam B(input.para()) ;
	cout << "Invarian mass of lepton pair " << M << endl ;
	cout << "Beam-up LO " << B.lo(500).u << endl;  
	cout << "Beam-down LO " << B.lo(500).d << endl;  
	cout << "Beam-up NLO " << B.nlo(500).u << endl;  
	cout << "Beam-down NLO " << B.nlo(500).d << endl;  

	// User code ends here
	//------------------------------------------------------- 
	t.stop();
	cout << "Time used : " << t.duration() << " secs" << endl;
	return 0;
}	
