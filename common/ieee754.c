#include <fpu_control.h>


#ifdef _FPU_SETCW
void __attribute__ ((constructor)) ieee754() 
          {
	    int cw;
	    cw=((_FPU_DEFAULT & ~_FPU_EXTENDED)|_FPU_DOUBLE) ;
            _FPU_SETCW (cw);
          }
#else

void __attribute__ ((constructor)) ieee754() 
          {
            (void) __setfpucw ((_FPU_DEFAULT & ~_FPU_EXTENDED)|_FPU_DOUBLE);
          }

#endif

