# This shell generates a grid for bottom production in p-barp collisions
# at the Tevatron. The grid is then convoluted with a non-perturbative
# fragmentation function and integrated over over y,
# so as to reproduce the B meson prediction published in PRL89 (2002) 122003.

# NB. In general, It is up to the user to make sure that the number and
# location of the points used for the grid produce an acceptable
# precision after interpolation and integration


#!/bin/sh

# output filename
file=Tevatron-1800

table=no

if [ $table = 'yes' ]; then 

# link to PDF files
../pdfdata/linkpdf.sh


# create the executable
(cd ../Linux ; makefonll )

# executable
FONLL=../Linux/fonll

# remove old files, if present
\rm $file.log $file.dat
touch $file.log  $file.dat

# pt values for the grid
for pt in 3. 5. 7. 10. 15. 20. 30. 40. 50. 75. 100. 125. 150. 200.
do

# y values for the grid
for ylab in 0 .25 .5 .75 1. 1.25
do

$FONLL << EOD >> $file.log
tmp$$
 1  900. 0 0 108 ! beam1: type, ener.,  nptype, ngroup, nset
 -1 900. 0 0 108 ! beam2: type, ener.,  nptype, ngroup, nset
  4.75 ! heavy quark mass
  -1       ! Lambda_5, <0 for default
  1 1 ! ren. and fact. scale factors
  $pt  $ylab ! pt,y_lab
 1        ! icalctype (1 ->FONLL, 2->massive)
EOD

cat tmp$$.out >> $file.dat
cat tmp$$.outlog >> $file.outlog
\rm tmp$$* 
\rm hdms-log.tmp  hdrsout.tmp phms-log.tmp hdres.tmp phres.tmp
done
done

\rm cor* cteq* *MRS* PION* vnv* LAC* GR* ft*
\rm $file.log $file.outlog

fi

# grid done. Interpolation and integration follows

##############################################################################
# This part of the shell convolutes the b-quark p_T distribution
# with a non-perturbative fragmentation function and integrates the
# grid. The reulting p_T distribution for the B mesons are in the file
# B-meson.top (the unconvoluted distribution is also included in the
# file)

# make the program interp
make

# pt distribution of B meson, integrated over y
file1=B-meson

./interp << EOD > /dev/null
$file1
2              ! icross:0->only interpol.,1->\int dydpt,2->\int dy,3->\int dpt
1               ! iwhatdpt: output as: 1 -> d/dpt, 2 -> d/dpt^2
5 30 26               ! ptmin,ptmax,ptpoints
0 1. 5            ! ymin,ymax,ypoints
1 0.375 29.1 0.    ! iwhatnp,br,par1,par2
4.75 900 900   ! qm,enh1,enh2
2               ! mult: moltiplicative factor for rapidity integrals
1 $file.dat  ! table format and filename
EOD

\rm $file1-interp.log    # the inputs above
\rm $file1-inputs.dat  # the input points and the convolution  
\rm $file1-fonllconv.dat # the interpolated points in the integration region
\rm $file1-fonll.dat # the interpolated points in the integration region

####################################################

exit
