#include <math.h>
#include <input.h>
#include <const.h>

using namespace constants;

ReadInput input;
// input SM parameters
namespace smPara{
	const int Nf = input.Nf;
	const double MW = input.MW;
	const double MZ = input.MZ;
	const double GF = input.GF;
}

// input Beam parameters
namespace beamPara{
	const double pTveto = input.pTveto;
	const double ECM = input.ECM;
	const double R = input.R;
}

// Derived SM parameters
namespace smPara{
	const double MW2=MW*MW; 
	const double MZ2=MZ*MZ;
	const double SW2=1.-MW2/MZ2;
	const double SW=sqrt(SW2);
	const double CW2=1.-SW2;	
	const double CW=sqrt(CW2);
	const double alpha = sqrt(2.)*MW2/pi*GF;
	const double alpha2 = alpha*alpha ;
}
