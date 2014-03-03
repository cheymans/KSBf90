import numpy as np


tiles_File = './tile_group.dat'

catalogue_Dir = '/home/cajd/Size_Measurement/PSF/Border_Removed/'
catalogue_File_Tail = '.rrg.shape.cat'

# Tile limits details which range of tile should be included #
tile_limits = [1, 80]

putHeader = True
Headers = (['Tile',
           'Number', 'Flux', 'Unknown', 'Mag',
           'Unknown','Unknown','Unknown','Unknown',
           'x', 'y', 'RA', 'Dec', 
           'Unknown', 'Unknown', 'Unknown', 'Unknown', 
           'Unknown', 'Unknown', 'Unknown', 'fr', 
           'FWHM', 'Unknown', 'nTot',
           'rg', 'e1_obs', 'e2_obs', 'e1_corr', 'e2_corr',
           'Pgamma', 'RRG Size (TrJ)', 'RRG Size (detJ)', 'KSB Size (TrJ)', 'KSB Size (detJ)'])
## Appended in this program [1]
##Start of original SExtractor Catalogue [2-5]
##[6-9]
#[10-13]
#[14-17]
#[18-21]
##End of SExtractor Catalogue [22-24]
##Start of KSBf90 addition

full_Catalogue_File = open(catalogue_Dir+'Full_Catalogue.cat', 'w')

if(putHeader):
    for header_index in range(0, len(Headers)):
        full_Catalogue_File.write('# '+str(header_index+1)+' '+Headers[header_index]+'\n')

totalSize = 0
for tile in range(tile_limits[0], tile_limits[1]+1):
    print 'Doing tile:', tile


    if(tile <10):
        stringTile = '0'+str(tile)
    else:
        stringTile = str(tile)
    # Read in Catalogue #
    Single_Cat_In = np.genfromtxt(catalogue_Dir+stringTile+catalogue_File_Tail)

    ### Add Extra Information ###
    Single_Cat = np.zeros((Single_Cat_In.shape[0], Single_Cat_In.shape[1]+1))
    Single_Cat[:,1:] = Single_Cat_In

    for i in range(0, Single_Cat.shape[0]):
        Single_Cat[i,0] = tile
    
    ### Selection Criteria must go here ###
    #######################################
    
    totalSize += Single_Cat.shape[0]

    np.savetxt(full_Catalogue_File, Single_Cat, fmt = '%10.4f')


full_Catalogue_File.close()

print 'Total of ', totalSize, ' items output'

print 'Finished without problems'

