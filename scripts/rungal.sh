#!/bin/sh

KSBDIR=/home/cajd/Size_Measurement/KSBf90/src/
DATADIR=/disk1/ps1/cech/STAGES/RRG_update/images/scifits/
STARCATDIR=/disk1/ps1/cech/STAGES/RRG_update/images/stars/starcats/
CATDIR=/disk1/ps1/cech/STAGES/RRG_update/images/SExcats/
OUTDIR=/home/cajd/Size_Measurement/PSF/
fitshead=a901-
fitstail=-sci.fits
starcattail=.chosen_stars.cat
sexcattail=-ab.cat

# first create image lists in 7 groups 
# Groups based on HST observation dates (see Fig 1 of Heymans et al 2008)

for obs_group in 1 2 3 4 5 6 7

do 
    # write out image location

    while read tile grp
    do
	if [ $grp = $obs_group ]; then

	    i=01

    # Run galcorrect on tile with correct obs_group PSF model
	    
	    $KSBDIR/gal_correct.a -image $DATADIR/$fitshead$tile$fitstail -in $CATDIR/$fitshead$tile$sexcattail -psf $OUTDIR/STAGES_PSF.$obs_group.dat -RRGpsf $OUTDIR/STAGES_RRG_PSF.$obs_group.dat -crit $STARCATDIR/$fitshead$i$starcattail -param KSBf90.param -out $OUTDIR/$tile.rrg.shape.cat -PSF_FWHM $OUTDIR/FWHM.$obs_group.dat

	fi
    done < tile_group.dat  
done

