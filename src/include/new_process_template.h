#ifndef NEW_PROCESS_H
#define NEW_PROCESS_H

#include <const.h>
#include <process.h>

using namespace smPara;

class new_process : public process{
	private :
	public :
		new_process() ;
		void calc_lo(), calc_nlo() ; //Calculate H0 and H1 
};
#endif
