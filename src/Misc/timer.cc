//  functions in timer class
//----------------------------------------------------------------
#include <misc.h>

void timer:: start() {start_time = clock() ; }
void timer:: stop() {stop_time = clock() ; }
double timer :: duration() {
	duration_time = stop_time - start_time ; 
	return double(duration_time)/((double)CLOCKS_PER_SEC) ;
}
