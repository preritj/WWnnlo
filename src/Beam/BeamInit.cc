#include <iostream> 
#include <cstdlib> 
#include <math.h>
#include <string>
#include <vector>
#include <algorithm>
#include <const.h>
#include <functions.h>
#include <Beam.h>

using namespace std;

Beam::Beam(){
	// Read input 
	vector<string> input_list ;
	input_list.clear() ;
	string input ;
	int comment_pos ;
	bool eof = false ;
	while(!eof) {
		getline(cin, input) ;
		comment_pos = input.find_first_of("#") ;
		if (comment_pos == 0) continue ;
		else if (comment_pos > 0) input.resize(comment_pos) ;
		if (input == "" || input[0] == ' ') eof = true ;
		else {
			input.erase(remove_if(input.begin(), input.end(), ::isspace), input.end());
			input_list.push_back(input) ;
		}
	}
	int i = 0;
	para.ECM = atof(input_list[i++].c_str()) ;
	char const* PDFname = input_list[i++].c_str();
	int const mem = (int)(atof(input_list[i++].c_str())+0.5) ;
	para.Nf = atof(input_list[i++].c_str()) ;
	para.pTveto = atof(input_list[i++].c_str()) ;
	para.mu = atof(input_list[i++].c_str()) ;
	para.R = atof(input_list[i].c_str()) ;

	// Initialize PDF	
	pdfini_(PDFname, &mem);
}

