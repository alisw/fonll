ifeq ($(shell uname),AIX)
F77=xlf -qextname -qflttrap=overflow:zerodivide:invalid:enable -O3 -qstrict
SYSOBJ=
AUTODBL=-qautodbl=dblpad
endif
ifeq ($(shell uname),SunOS)
F77= f77 -fnonstd
SYSOBJ=
endif
ifeq ($(shell uname),Linux)
F77= g77 -Wall -fno-automatic -ffortran-bounds-check
DEBUG=-ggdb
#OPT=-O
SYSOBJ=trapfpe.o
endif
ifeq ($(shell uname),Darwin)
F77= g77 -Wall -fno-automatic -ffortran-bounds-check
OPT=-O
#SYSOBJ=trapfpe.o
endif
ifeq ($(shell uname),HP-UX)
F77= g77 -Wall
SYSOBJ=
endif

FF=$(F77) $(OPT) $(DEBUG) $(AUTODBL) 


INTPACK = fakenag.o vegas_p.o
STRF =  fonstru_mlm.o fonllpdf.o
STRFMLM = hvqpdfpho.o dummies.o elpdf_user.o
STRFPDL = jetpdflib.o elpdf_user.o dummies.o
#LIBSPDL = /afs/cern.ch/asis/i386_redhat72/cern/pro/lib/libpdflib804.a \
#          /afs/cern.ch/asis/i386_redhat72/cern/pro/lib/libmathlib.a \
#          /afs/cern.ch/asis/i386_redhat72/cern/pro/lib/libkernlib.a
LIBSPDL = -L/cern/pro/lib -lpdflib -lpacklib -lmathlib -lkernlib
COMMON= cernlib.o dfint.o gauleg.o fww_ww.o mlmident.o locatd.o dcspln.o\
          strf.o getstr.o spmat1.o fonfra_p.o istrl.o toend.o
#LIBS    = -L/cern/pro/lib -lpacklib -lmathlib -lkernlib
INCLUDE = xsect.h

VPATH=./:../hdmassive:../hdresummed:../phmassive:../phresummed:../common

alfas.o fbm.o fonfra_p.o fonstru.o fonstru_mlm.o \
                        input.o main_p.o testtstru.o: xsect.h

%.o: %.f
	$(FF) -c $<
%.o: %.c
	$(CC) $(DEBUG) -c $^

DRIVER=fonll0.o
DRIVER1=fonll.o
DRIVER2=fonlltest.o
DRIVER3=fonllgrid.o

HDM=hdms.o hvqcrf.o rozlimqa.o rozlimgg.o rozlimqg.o

HDR=main_p.o input.o matrix_p.o stru.o fbm.o

PHM=phms.o phocrs.o rozlimpg.o rozlimpc.o  rozlimpn.o  rozlimpa.o

PHR=photon_wps8.o fbm-gamma.o

fonll : $(DRIVER1) $(DRIVER) $(HDM) $(HDR) $(PHM) $(PHR) $(STRF) $(STRFMLM) \
        $(COMMON) $(INTPACK) $(SYSOBJ)
	$(FF) $^ $(LIBS) -o $@

fonlltest : $(DRIVER2) $(DRIVER) $(HDM) $(HDR) $(PHM) $(PHR) $(STRF) \
        $(STRFMLM) $(COMMON) $(INTPACK) $(SYSOBJ)
	$(FF) $^ $(LIBS) -o $@

fonllgrid: $(DRIVER3) pttoz.o $(DRIVER) $(HDM) $(HDR) $(PHM) $(PHR) $(STRF) \
        $(STRFMLM) $(COMMON) $(INTPACK) $(SYSOBJ)
	$(FF) $^ $(LIBS) -o $@

fonllpdl : $(DRIVER1) $(DRIVER) $(HDM) $(HDR) $(PHM) $(PHR) $(STRF) \
        $(STRFPDL) $(COMMON) $(INTPACK) $(SYSOBJ)
	$(FF) $^ $(LIBSPDL) -o $@

fragmfonll: fragmfonll.o bispline.o pttoz.o vegas_p.o toend.o cernlib.o
	$(FF) $^ -o $@

testpdl: testpdl.o
	$(FF) $^ $(LIBSPDL) -o $@

testmlm: testmlm.o $(STRFMLM)
	$(FF) $^ $(LIB) -o $@

tgz:
	tar zcvf massive.tgz ./ --exclude "*core" \
         --exclude massive.tgz \
         --exclude "*SunOS*" --exclude "*AIX*" --exclude "*Linux*" \
         --exclude "*.o" --exclude "*~" --exclude "*#" --exclude "*xsectmlm" --
clean:
	-rm core *.o *~



