### Python code to remove overlap regions according to the 0/1 value in the border files ###
import pyfits as pf
import numpy as np

import pylab as pl

border_Dir = '/disk1/ps1/cech/STAGES/RRG_update/images/border/'
border_Head = 'a901-'
border_Tail = '-border.fits.gz'

catalogue_Dir = '/home/cajd/Size_Measurement/PSF/'
catalogue_Tail = '.rrg.shape.cat'
x_pos_index = 9
y_pos_index = 10

output_Dir = '/home/cajd/Size_Measurement/PSF/Border_Removed/'

tile_Range = [1, 80]


x_pos_index -= 1; y_pos_index -= 1 ##Convert to python zero based##
for tile in range(tile_Range[0], tile_Range[1]+1):
    if(tile < 10):
        border_File = border_Head+'0'+str(tile)+border_Tail
        catalogue_File = '0'+str(tile)+catalogue_Tail
    else:
        border_File = border_Head+str(tile)+border_Tail
        catalogue_File = str(tile)+catalogue_Tail
    

    ##Read in Catalogue##
    Cat = np.genfromtxt(catalogue_Dir+catalogue_File)
    Discard_Flag = np.zeros(Cat.shape[0]) #0 = Keep, 1 = Discard

    ##Read in FITS border file##
    border_Data = pf.getdata(border_Dir+border_File)

    print '!-------------------------------------------------------------------'
    print border_Data.shape

    print 'Filename:', catalogue_File
    print 'Pos 1:', Cat[1,x_pos_index], Cat[1,y_pos_index]

    print 'Max/Min x:', np.amax(Cat[:,x_pos_index]), np.amin(Cat[:,x_pos_index])
    print 'Max/Min y:', np.amax(Cat[:,y_pos_index]), np.amin(Cat[:,y_pos_index])
    print '!-------------------------------------------------------------------'

    for gal in range(0, Discard_Flag.shape[0]):
        if(border_Data[np.round(Cat[gal,y_pos_index]).astype(int)-1,np.round(Cat[gal,x_pos_index]).astype(int)-1] == 1): #Change 1 here to flag for masked area#
            ## Discard ##
            Discard_Flag[gal] = 1

    number_To_Keep = Discard_Flag.shape[0]-np.sum(Discard_Flag)
    print ' Keeping ', number_To_Keep, ' of', Discard_Flag.shape[0], 'galaxies'

    ## Construct New Catalogue ###
    new_Catalogue = np.zeros((number_To_Keep, Cat.shape[1]))
    counter = 0
    for gal in range(0, Discard_Flag.shape[0]):
        if(Discard_Flag[gal] == 0):
            new_Catalogue[counter,:] = np.copy(Cat[gal,:])
            counter += 1

    print 'Counter:', counter, number_To_Keep
    ### Output ###
    np.savetxt(output_Dir+catalogue_File, new_Catalogue, fmt = '%10.4f')

    ### Plot ###
    '''
    f = pl.figure()
    ax = f.add_subplot(1,1,1)

    ax.imshow(border_Data, interpolation = None)
    for gal in range(0, new_Catalogue.shape[0]):
        print 'Plotting:', gal
        ax.plot(np.round(new_Catalogue[gal,x_pos_index]).astype(int)-1, np.round(new_Catalogue[gal,y_pos_index]).astype(int)-1, marker='o', markerfacecolor='red')

    pl.show()
    '''
