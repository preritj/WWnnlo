#include <iostream> 
#include <cstdlib> 
#include <math.h>
#include <string>
#include <vector>
#include <algorithm>
#include <const.h>
#include <functions.h>
#include <Beam.h>
#include <input.h>

using namespace std;

Beam::Beam(InputPara input){
	para.ECM = input.ECM ;
	const char* PDFname = input.PDFname;
	int  mem = input.mem ;
	para.Nf = input.Nf ;
	para.pTveto = input.pTveto ;
	para.R = input.R;
	para.mu= input.pTveto; // default scale = pT veto
	// Initialize PDF	
	pdfini_(PDFname, &mem);
}

