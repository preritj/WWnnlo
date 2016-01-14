#ifndef FUNCTIONS_H
#define FUNCTIONS_H

// PDF functions 
void evolvePDF(double , double , double* );
double alpha_s (double);

// Hard functions
double H1(double, double) ;

//Linking to PDFs (Fortran code for lhapdf) 
extern "C" 
{ 
 void pdfini_(char const*, int const*); 
 void evolvepdf_(double&, double&, double*  );
 double alphaspdf_(double&) ;
}

#endif
