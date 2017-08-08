# This shell generates a grid for charm electroproduction at HERA using
# the ETAG33 cuts of the H1 experiments. The grid is then convoluted
# with a Peterson non-perturbative fragmentation function and integrated
# over pt or y to produce results to be compared with the published data.
# It is up to the user to make sure that the number and location of the 
# points used for the grid produce an acceptable precision after
# interpolation and integration


#!/bin/sh

# output filename
file=H1-etag33

# this switch allows one to re-use the same table (set table=no) for different
# integrations/ convolutions, without regenerating it (this is the
# time-consuming step)
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
for pt in 2. 2.5 3. 3.5 4. 4.5 5. 5.5 6. 7. 8. 9. 10. 11. 12. 13. 14. 15. 20
do

# y values for the grid
for ylab in -1.8 -1.7 -1.6 -1.55 -1.5 -1.45 -1.4 -1.35 -1.3 -1.2 -1.1 -1 -.75 -.5 -.25 0 .25 .5 .75 1 1.25 1.5 1.75 2.
do

$FONLL << EOD >> $file.log
tmp$$
 1  820. 0 0 101 ! beam1: type, ener.,  nptype, ngroup, nset
 5  27.5 0 0 42 ! beam2: type, ener.,  nptype, ngroup, nset
  1.5 ! heavy quark mass
  -1       ! Lambda_5, <0 for default
  1 1 ! ren. and fact. scale factors
  $pt  $ylab ! pt,y_lab
 1        ! icalctype (1 ->FONLL, 2->massive)
1	! 1-> cut on photon virtuality
.1 ! upper scale in WW (GeV)
.29 .62 ! zmin, zmax
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

# make the program interp
make

# y distribution, integrated over pt 
file1=$file-y

./interp << EOD > /dev/null
$file1
3               ! icross:0->only interpol.,1->\int dydpt,2->\int dy,3->\int dpt
1               ! iwhatdpt: output as: 1 -> d/dpt, 2 -> d/dpt^2
5 10.5 11               ! ptmin,ptmax,ptpoints
-1.5 1.5 31             ! ymin,ymax,ypoints
2 1. 0.02 0.    ! iwhatnp,br,par1,par2
1.5 820. 27.5   ! qm,enh1,enh2
1               ! mult: moltiplicative factor for rapidity integrals
1 $file.dat  ! table format and filename
EOD

\rm $file1-interp.log    # the inputs above
\rm $file1-inputs.dat  # the input points and the convolution  
\rm $file1-fonllconv.dat # the interpolated points in the integration region
\rm $file1-fonll.dat # the interpolated points in the integration region

####################################################

# pt distribution, integrated over y
file1=$file-pt

./interp << EOD > /dev/null
$file1
2               ! icross:0->only interpol.,1->\int dydpt,2->\int dy,3->\int dpt
1               ! iwhatdpt: output as: 1 -> d/dpt, 2 -> d/dpt^2
2 12 11               ! ptmin,ptmax,ptpoints
-1.5 1.5 31             ! ymin,ymax,ypoints
2 1. 0.02 0.    ! iwhatnp,br,par1,par2
1.5 820. 27.5   ! qm,enh1,enh2
1               ! mult: moltiplicative factor for rapidity integrals
1 $file.dat  ! table format and filename
EOD

\rm $file1-interp.log    # the inputs above
\rm $file1-inputs.dat  # the input points and the convolution  
\rm $file1-fonllconv.dat # the interpolated points in the integration region
\rm $file1-fonll.dat # the interpolated points in the integration region

# the results in the files $file1-y.top and $file1-pt.top should be
# multiplied by 470 (pb->nb, branching ratio, c+cbar) in order to be
# compared to those in Fig. 8 of Frixione-Nason

exit
