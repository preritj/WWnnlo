#include <iostream> 
#include <cstdlib> 
#include <math.h>
#include <string>
#include <vector>
#include <algorithm>
#include <const.h>
#include <PDF.h>
#include <Beam.h>
#include <const.h>
#include <input.h>

using namespace std;

int Nf;
double MW, MZ, GF ;
double pTveto, R, ECM;

ReadInput::ReadInput(){
	// Read input 
	vector<string> input_list ;
	input_list.clear() ;
	string inputline ;
	int comment_pos ;
	bool eof = false ;
	while(!eof) {
		getline(cin, inputline) ;
		comment_pos = inputline.find_first_of("#") ;
		if (comment_pos == 0) continue ;
		else if (comment_pos > 0) inputline.resize(comment_pos) ;
		if (inputline == "" || inputline[0] == ' ') eof = true ;
		else {
			inputline.erase(remove_if(inputline.begin(), inputline.end(), ::isspace),
							 inputline.end());
			input_list.push_back(inputline) ;
		}
	}
	int i = 0;
	ECM = atof(input_list[i++].c_str()) ;
	const char* PDFname = input_list[i++].c_str();
	const int mem = (int)(atof(input_list[i++].c_str())+0.5) ;
	Nf = atof(input_list[i++].c_str()) ;
	pTveto = atof(input_list[i++].c_str()) ;
	R = atof(input_list[i++].c_str()) ;
	MW = atof(input_list[i++].c_str()) ;
	MZ = atof(input_list[i++].c_str()) ;
	GF = atof(input_list[i].c_str()) ;
	SMinit();
	// Initialize PDF	
	pdfini_(PDFname, &mem);
}



