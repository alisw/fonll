ifeq ($(shell uname),Linux)
FC=g77 -g -Wall -fno-automatic
FFLAGS=
LIBS = -L/cern/pro/lib -lpacklib -lmathlib -lkernlib \
       
endif

interp : interp.o bispline.o d01fcf.o vegas_p.o Makefile
	$(FC) interp.o bispline.o d01fcf.o vegas_p.o $(LIBS) -o $@

.f.o:
	$(FC) -c $*.f
