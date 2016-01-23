#include <math.h>
#include <input.h>
#include <const.h>

double MW2, MZ2;
double SW, CW, SW2, CW2;
double alpha, alpha2 ;

void ReadInput::SMinit(){
	MW2=MW*MW; 
	MZ2=MZ*MZ;
	SW2=1.-MW2/MZ2;
	SW=sqrt(SW2);
	CW2=1.-SW2;	
	CW=sqrt(CW2);
	alpha = sqrt(2.)*MW2/pi*GF;
	alpha2 = alpha*alpha ;
}
