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

Beam::Beam(){
	para.mu= pTveto; // default scale = pT veto
}

