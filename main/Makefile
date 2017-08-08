ifeq ($(shell uname),Linux)
#FC= g77 -Wall -fno-automatic -ffortran-bounds-check
# Make sure that, with gfortran, -O is used together together with -fno-automatic 
# (-O gives about 30% improvement in time, but is broken without -fno-automatic
FC= gfortran -Wall -fno-automatic
#FF=/usr/bin/g95 -Wall -fno-automatic
#FF=/usr/local/ifort/bin/ia32/ifort -save
#DEBUG=-ggdb
OPT=-O
#STATIC=-static
SYSOBJ=trapfpe.o
endif

ifeq ($(shell uname),Darwin)
#FC= g77 -Wall -fno-automatic -ffortran-bounds-check
FC= gfortran -Wall -fno-automatic
OPT=-O
STATIC=
endif

FF=$(FC) $(OPT) $(DEBUG) $(AUTODBL) 

INTPACK = fakenag.o vegas_p.o
STRF =  fonstru_mlm.o fonllpdf.o eks98.o eps09.o
STRFMLM = hvqpdfpho.o dummies.o elpdf_user.o
STRFPDL = jetpdflib.o elpdf_user.o dummies.o
STRFLHA = lhapdfif.o

LIBSPDL = -L/cern/pro/lib -lpdflib -lpacklib -lmathlib -lkernlib
# the following line needs lhapdf-config to be correctly in PATH
LIBSLHA = $(STATIC) -L$(shell lhapdf-config --prefix)/lib -lLHAPDF
COMMON= cernlib.o dfint.o gauleg.o fww_ww.o mlmident.o \
          strf.o getstr.o spmat1.o fonfra_p.o istrl.o toend.o dcspln.o locatd.o
INCLUDE = xsect.h

VPATH=./:../hdmassive:../hdresummed:../phmassive:../phresummed:../common

alfas.o fbm.o fonfra_p.o fonstru.o fonstru_mlm.o \
                        input.o main_p.o testtstru.o: xsect.h

%.o: %.f
	$(FF) -c $<
%.o: %.c
	$(CC) $(DEBUG) -c $^

DRIVER=fonll0.o
DRIVER1= fonll.o
DRIVER2= fonlltest.o

HDM=hdms.o hvqcrf.o rozlimqa.o rozlimgg.o rozlimqg.o

HDR=main_p.o input.o matrix_p.o stru.o fbm.o struv_plus_avgo.o

PHM=phms.o phocrs.o rozlimpg.o rozlimpc.o  rozlimpn.o  rozlimpa.o

PHR=photon_wps8.o fbm-gamma.o

fonll : $(DRIVER1) $(DRIVER) $(HDM) $(HDR) $(PHM) $(PHR) $(STRF) $(STRFMLM) \
        $(COMMON) $(INTPACK) $(SYSOBJ)
	$(FF) $^ $(LIBS) -o $@

fonlltest : $(DRIVER2) $(DRIVER) $(HDM) $(HDR) $(PHM) $(PHR) $(STRF) \
        $(STRFMLM) $(COMMON) $(INTPACK) $(SYSOBJ)
	$(FF) $^ $(LIBS) -o $@

fonllpdl : $(DRIVER1) $(DRIVER) $(HDM) $(HDR) $(PHM) $(PHR) $(STRF) \
        $(STRFPDL) $(COMMON) $(INTPACK) $(SYSOBJ)
	$(FF) $^ $(LIBSPDL) -o $@

fonlllha : $(DRIVER1) $(DRIVER) $(HDM) $(HDR) $(PHM) $(PHR) $(STRF) \
        $(STRFLHA) $(COMMON) $(INTPACK) $(SYSOBJ)
	$(FF) $^ $(LIBSLHA) -o $@

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



