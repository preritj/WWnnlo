//  timer class
#include <time.h> 
class timer
{
	private :
		clock_t start_time, stop_time, duration_time ;
	public : 
		void start(), stop() ;
		double duration() ;
};


// function for Integration 
void NIntegrate(double, double, void*, 
				double (*integrand)(double, void*), 
				double& , double& , int ) ;

