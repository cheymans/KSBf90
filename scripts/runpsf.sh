#!/bin/sh

KSBDIR=/home/cech/KSBf90
DATADIR=/disk1/ps1/cech/STAGES/RRG_update/images/scifits/
STARCATDIR=/disk1/ps1/cech/STAGES/RRG_update/images/stars/starcats/
CATDIR=/disk1/ps1/cech/STAGES/RRG_update/images/SExcats/
OUTDIR=/disk1/ps1/cech/STAGES/RRG_update/OUTDATA/PSF
fitshead=a901-
fitstail=-sci.fits
starcattail=.chosen_stars.cat
sexcattail=-ab.cat

# first create image lists in 7 groups 
# Groups based on HST observation dates (see Fig 1 of Heymans et al 2008)

for obs_group in 1 #2 3 4 5 6 7
do 

    # write out image location

    while read tile grp
    do
	if [ $grp = $obs_group ]; then

	    echo $DATADIR/$fitshead$tile$fitstail
	fi
    done < tile_group.dat  >  IM_CAT_LISTS/image.list.$obs_group

    # write out catalogue location

    while read tile grp
    do
	if [ $grp = $obs_group ]; then

	    echo $CATDIR/$fitshead$tile$sexcattail

	fi
    done < tile_group.dat  >  IM_CAT_LISTS/cat.list.$obs_group


    # count the number of images in each group
    COUNTER=`wc -l IM_CAT_LISTS/cat.list.$obs_group | awk '{print $1}'`
    
    i=01

    # RUN psffit for each group to create PSF models

    $KSBDIR/psffit.a -nimage $COUNTER -image IM_CAT_LISTS/image.list.$obs_group -in IM_CAT_LISTS/cat.list.$obs_group -crit $STARCATDIR/$fitshead$i$starcattail -param KSBf90.param -pgopen /xwin -out $OUTDIR/STAGES_PSF.$obs_group.dat

done

