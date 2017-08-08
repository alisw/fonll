#!/bin/bash

#       	 hvq       lha
#	cteq6m	 131	10050
#	AFG MC	  42	362
#	SMRS-P 2  32	232


# Chose the desired PDF library by commenting/uncommenting
# "hvq" uses hvqpdfpho.f, "lha" uses the external LHAPDF library
#pdflib='hvq'
pdflib='lha'

if [[ $pdflib == 'hvq' ]]; then
  proton=131
  photon=42
  pion=32
  outfile=testfonll2-hvq
  FONLL=../`uname`/fonll
  ../pdfdata/linkpdf.sh
else
  proton=10050
  photon=363
  pion=232
  outfile=testfonll2-lha
  FONLL=../`uname`/fonlllha
fi

\rm $outfile.outlog
echo $FONLL


# electroproduction at HERA

# for ffact in 0.5 1 2
for ffact in 2
do
# for fren in 0.5 1 2
for fren in 0.5
do
# for pt in 0.5 5 20
for pt in 20
do
# for y in -1 0 1
for y in 1
do
$FONLL <<xxx
$outfile
 1  920. 0 0 $proton ! beam1: type, ener.,  nptype, ngroup, nset
 5  27.5 0 0 $photon ! beam2: type, ener.,  nptype, ngroup, nset
  1.5 ! heavy quark mass
 -1    ! Lambda_5, <=0 for default
  $ffact $fren !  fact. and ren. scale factors
  $pt $y ! pt,y_lab
 1   ! icalctype
 1    ! itype ww (1-4)
  1.   ! effective ww scale
  0.2  0.8 ! zminww,zmaxww
xxx
done
done
done
done


# same thing, exchanging the beams, pt=20, fact e ren=2 e 0.5, y=1

$FONLL <<xxx
$outfile
 5  27.5 0 0 $photon ! beam1: type, ener.,  nptype, ngroup, nset
 1  920. 0 0 $proton ! beam2: type, ener.,  nptype, ngroup, nset
  1.5 ! heavy quark mass
 -1    ! Lambda_5, <=0 for default
  2 0.5  !  fact. and ren. scale factors
  20 -1 ! pt,y_lab
 1   ! icalctype
 1    ! itype ww (1-4)
  1.   ! effective ww scale
  0.2  0.8 ! zminww,zmaxww
xxx


# photon at fixed target

# for pt in 0.5 5 10
for pt in 5
do
# for y in 1 2 3
for y in 2
do
$FONLL <<xxx
$outfile
 4  200. 0 0 $photon ! beam1: type, ener.,  nptype, ngroup, nset
 0  0    0 0 $proton ! beam2: type, ener.,  nptype, ngroup, nset
  1.5 ! heavy quark mass
 -1    ! Lambda_5, <=0 for default
  0.5  2 !  fact. and ren. scale factors
  $pt $y ! pt,y_lab
 1   ! icalctype
xxx
done
done

# bottom at Tevatron (900+900):

# for pt in 2 5 10 20
for pt in 20
do
# for y in -1 0 1
for y in 1
do
$FONLL <<xxx
$outfile
 1  900. 0 0 $proton ! beam1: type, ener.,  nptype, ngroup, nset
 -1 900. 0 0 $proton ! beam2: type, ener.,  nptype, ngroup, nset
  5    ! heavy quark mass
 -1    ! Lambda_5, <=0 for default
  2.  2. !  fact. and ren. scale factors
  $pt $y ! pt,y_lab
 1   ! icalctype
xxx
done
done

# bottom at Tevatron (980+980):

# for pt in 2 5 10 20
for pt in 20
do
# for y in -1 0 1
for y in 1
do
$FONLL <<xxx
$outfile
 1  980. 0 0 $proton ! beam1: type, ener.,  nptype, ngroup, nset
 -1 980. 0 0 $proton ! beam2: type, ener.,  nptype, ngroup, nset
  5    ! heavy quark mass
 -1    ! Lambda_5, <=0 for default
  2.  2. !  fact. and ren. scale factors
  $pt $y ! pt,y_lab
 1   ! icalctype
xxx
done
done

# charm at Tevatron (900+900):

# for pt in 2 5 10 20
for pt in 10
do
# for y in -1 0 1
for y in 1
do
$FONLL <<xxx
$outfile
 1  900. 0 0 $proton ! beam1: type, ener.,  nptype, ngroup, nset
 -1 900. 0 0 $proton ! beam2: type, ener.,  nptype, ngroup, nset
  1.5    ! heavy quark mass
 -1    ! Lambda_5, <=0 for default
  2.  1. !  fact. and ren. scale factors
  $pt $y ! pt,y_lab
 1   ! icalctype
xxx
done
done


# bottom at LHC (7000+7000):

# for pt in 2 5 10 20
for pt in 5
do
# for y in -1 0 1
for y in 2
do
$FONLL <<xxx
$outfile
 1  7000. 0 0 $proton ! beam1: type, ener.,  nptype, ngroup, nset
 1 7000. 0 0 $proton ! beam2: type, ener.,  nptype, ngroup, nset
  5    ! heavy quark mass
 -1    ! Lambda_5, <=0 for default
  2.  2. !  fact. and ren. scale factors
  $pt $y ! pt,y_lab
 1   ! icalctype
xxx
done
done

# charm at LHC (7000+7000):

# for pt in 2 5 10 20
for pt in 5
do
# for y in -1 0 1
for y in 2
do
$FONLL <<xxx
$outfile
 1  7000. 0 0 $proton ! beam1: type, ener.,  nptype, ngroup, nset
 1 7000. 0 0 $proton ! beam2: type, ener.,  nptype, ngroup, nset
  1.5    ! heavy quark mass
 -1    ! Lambda_5, <=0 for default
  2.  1. !  fact. and ren. scale factors
  $pt $y ! pt,y_lab
 1   ! icalctype
xxx
done
done


# pion at fixed target

# for pt in 0.5 1 3 10
for pt in 3
do
# for y in 1 2 3 4
for y in 2
do
$FONLL <<xxx
$outfile
3  800 0 0 $pion   ! 2,2,3 in pdflib
0  0   0 0 $proton
1.5    ! heavy quark mass
-1    ! Lambda_5, <=0 for default
1.5 2   !  fact. and ren. scale factors
$pt $y   ! pt,y_lab
1     ! icalctype
xxx
done
done

\rm hdms-log.tmp  hdrsout.tmp phms-log.tmp hdres.tmp phres.tmp
\rm *.tmp $outfile'.out' $outfile'fonll.log'
\rm cor* ctq* cteq* *MRS* PION* vnv* LAC* GR* ft* alf*

exit
