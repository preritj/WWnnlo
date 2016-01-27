// QNG (quadrature) method is used if Ndim is not specified
#include <iostream>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
using namespace std;

void NIntegrate(double lower_limit, double upper_limit, 
		 void *params_ptr, double (*integrand)
		(double, void*), double& result, 
		double& error, int verbose )
{
	gsl_error_handler_t * old_handler=gsl_set_error_handler_off();
	gsl_integration_workspace *w 
		= gsl_integration_workspace_alloc (50000);

	double abs_error = 5.0e-3;	/* to avoid round-off problems */
	double rel_error = 1.0e-3;	/* the result will usually be much better */
	size_t Neval ;

	gsl_function My_function;

	My_function.function = integrand;
	My_function.params = params_ptr;
	
	int status=1;
    while(status) {
        status= gsl_integration_qags (&My_function, lower_limit, 
	 	upper_limit, abs_error, rel_error, 10000, w, &result, &error);
	 	if (status) rel_error*=2.0 ;
	 	if (verbose == 1) {
	 		printf ("result = % .6f error = % .6f \n", result, error);
	 		if (status) cout << "Increasing tolerance " << endl ;
	 		} ;
	 	}
    gsl_integration_workspace_free (w);
    gsl_set_error_handler(old_handler);
} 
