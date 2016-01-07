This program requires following programs to be pre-installed : 
* gcc
* gfortran
* gsl
* LHApdf (fortran version) 

----------------------------------------------------------------- 


In particular, $PATH variable should include 'bin' directories of gsl and lhapdf
while $LD_LIBRARY_PATH variable should include 'lib' directories of gsl and lhapdf

This can be implemented by adding following lines in ~/.bashrc : 

export PATH = $HOME/Apps/gsl/Install/bin : $HOME/Apps/LHAPDF/Install/bin : $PATH 

export LD_LIBRARY_PATH = $HOME/Apps/gsl/Install/lib : $HOME/Apps/LHAPDF/Install/lib : $LD_LIBRARY_PATH 

followed by 

source ~/.bashrc 

----------------------------------------------------------------- 

To compile, do : 

  make 

----------------------------------------------------------------- 


To recompile everything from scratch, do : 

  make clean 
  
  make 



The only files you need to edit can be found in src/User 

----------------------------------------------------------------- 

To run the program after compilation, go to directory bin 

Edit input file 'infile', and then do : 

./Beam < infile

