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
	input.ECM = atof(input_list[i++].c_str()) ;
	input.PDFname = input_list[i++].c_str();
	input.mem = (int)(atof(input_list[i++].c_str())+0.5) ;
	input.Nf = atof(input_list[i++].c_str()) ;
	input.pTveto = atof(input_list[i++].c_str()) ;
	input.R = atof(input_list[i].c_str()) ;
	input.MW = atof(input_list[i].c_str()) ;
	input.MZ = atof(input_list[i].c_str()) ;
	input.GF = atof(input_list[i].c_str()) ;
}

