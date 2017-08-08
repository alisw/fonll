#!/bin/bash

cat <<EOF > tmp.f
      call PRNTSF
      end
      subroutine elpdf_user
      end
EOF

g77 ../common/hvqpdfpho.f ../common/dummies.f tmp.f

cat <<EOF > hvqpdfpho.txt
Numbering scheme for the Parton Distribution Function sets contained
in the hvqpdfpho.f library.

NB. When doing electroproduction the user must chose a PHOTON PDF set;
    the FONLL program takes care of computing the corresponding electron
    parton density.
 
EOF
 


./a.out >> hvqpdfpho.txt

\rm a.out tmp.f
