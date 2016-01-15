#ifndef INPUT_H
#define INPUT_H

struct InputPara{
	double ECM;		     // CoM energy
	double pTveto;   	 // pT-veto
	double R; 			 // jet radius
	int Nf; 			 // number of flavors 
	const char* PDFname; // PDF name
	int mem ;		 // PDF member
	double MW, MZ, GF;	 // EW input parameters 
};


class ReadInput{
	private :
		InputPara input;	
	public :
		ReadInput();
		InputPara para() {return input;}; 
};

#endif