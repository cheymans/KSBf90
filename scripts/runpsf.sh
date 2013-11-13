#!/bin/sh

KSBDIR=/home/cajd/Size_Measurements/src/
DATADIR=/disk2/cech/STAGES_image_analysis/images/scifits/
STARCATDIR=/disk2/cech/STAGES_image_analysis/stars/starcats/
CATDIR=/disk2/cech/STAGES_image_analysis/images/SExcats/
fitshead=a901-
fitstail=-sci.fits
starcattail=.chosen_stars.cat
sexcattail=-ab.cat

# first create image lists

#for file in $DATADIR/*
#do
#echo $file 
#done > image.list

#for file in $CATDIR/*$sexcattail
#do 
#echo $file
#done > cat.list

i=01
#$KSBDIR/psffit.a -nimage 1 -image $DATADIR/$fitshead$i$fitstail -in $CATDIR/$fitshead$i$sexcattail -crit $STARCATDIR/$fitshead$i$starcattail -param KSBf90.param -pgopen /xwin -out PSF.dat

#$KSBDIR/psffit.a -nimage 3 -image image.list -in cat.list -crit $STARCATDIR/$fitshead$i$starcattail -param KSBf90.param -pgopen /xwin -out PSF.dat

# single image - hardwired to disk1 until disk2 faults corrected

$KSBDIR/psffit.a -nimage 1 -image /home/cech/a901-01-sci.fits -in /home/cajd/Size_Measurements/test_data/a901-01-ab.cat -crit  /home/cajd/Size_Measurements/test_data/a901-01.chosen_stars.cat -param KSBf90.param -pgopen /xwin -out PSF.dat