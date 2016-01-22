#include <iostream> 
#include <cstdlib> 
#include <math.h>
#include <string>
#include <vector>
#include <algorithm>
#include <const.h>
#include <PDF.h>
#include <Beam.h>
#include <input.h>

using namespace std;

Beam::Beam(InputPara input){
	para.ECM = input.ECM ;
	para.Nf = input.Nf ;
	para.pTveto = input.pTveto ;
	para.R = input.R;
	para.mu= input.pTveto; // default scale = pT veto
}

