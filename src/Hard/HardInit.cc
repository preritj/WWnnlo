#include <iostream> 
#include <cstdlib> 
#include <math.h>
#include <const.h>
#include <input.h>
#include <Hard.h>

using namespace std;

Hard::Hard(InputPara input){
	MW = input.MW;
	MZ = input.MZ;
	GF = input.GF;
	MW2=MW*MW; MZ2=MZ*MZ;
	SW2=1.-MW2/MZ2, SW=sqrt(SW2);
	CW2=1.-SW2; CW=sqrt(CW2);
}

