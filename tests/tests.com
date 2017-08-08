#!/bin/bash

FONLL=../`uname`/fonll

../pdfdata/linkpdf.sh

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
testfonll
 1  920. 0 0 131 ! beam1: type, ener.,  nptype, ngroup, nset
 5  27.5 0 0 42 ! beam2: type, ener.,  nptype, ngroup, nset
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
testfonll
 5  27.5 0 0 42 ! beam1: type, ener.,  nptype, ngroup, nset
 1  920. 0 0 131 ! beam2: type, ener.,  nptype, ngroup, nset
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
testfonll
 4  200. 0 0 42 ! beam1: type, ener.,  nptype, ngroup, nset
 0  0    0 0 131 ! beam2: type, ener.,  nptype, ngroup, nset
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
testfonll
 1  900. 0 0 131 ! beam1: type, ener.,  nptype, ngroup, nset
 -1 900. 0 0 131 ! beam2: type, ener.,  nptype, ngroup, nset
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
testfonll
 1  980. 0 0 131 ! beam1: type, ener.,  nptype, ngroup, nset
 -1 980. 0 0 131 ! beam2: type, ener.,  nptype, ngroup, nset
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
testfonll
 1  900. 0 0 131 ! beam1: type, ener.,  nptype, ngroup, nset
 -1 900. 0 0 131 ! beam2: type, ener.,  nptype, ngroup, nset
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
testfonll
3  800 0 0 32   ! 2,2,3 in pdflib
0  0   0 0 131
1.5    ! heavy quark mass
-1    ! Lambda_5, <=0 for default
1.5 2   !  fact. and ren. scale factors
$pt $y   ! pt,y_lab
1     ! icalctype
xxx
done
done

\rm hdms-log.tmp  hdrsout.tmp phms-log.tmp hdres.tmp phres.tmp
\rm testfonll*.tmp testfonll.out testfonllfonll.log
\rm cor* ctq* cteq* *MRS* PION* vnv* LAC* GR* ft* alf*

exit
