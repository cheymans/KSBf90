import pylab as pl
import numpy as np


file_Ellip = '../Corrected_Ellipticity_Comparison.dat'
file_Size = '../Corrected_Size_Comparison.dat'

E_Limits = [-0.4, 0.4]

# Column Numbers for Comparison #
RRG_E_Col = [0,1]
KSB_E_Col = [2,3]
RRG_S_Col = [0,1]
KSB_S_Col = [2,3]

yeqx = np.linspace(E_Limits[0], E_Limits[1], 50)

###Produce Ellipticity Plot##

#Read in#
E_Input = np.genfromtxt(file_Ellip)
RRG_E = np.zeros((E_Input.shape[0],2))
KSB_E = np.zeros((E_Input.shape[0],2))

RRG_E[:,0] = E_Input[:,RRG_E_Col[0]]; RRG_E[:,1] = E_Input[:,RRG_E_Col[1]]
KSB_E[:,0] = E_Input[:,KSB_E_Col[0]]; KSB_E[:,1] = E_Input[:,KSB_E_Col[1]]

fE = pl.figure(0, (6,2*6))
ax = fE.add_subplot(2,1,1)

ax.set_xlabel('RRG Ellipticity')
ax.set_ylabel('KSB Ellipticity')

ax.set_title('e(1)')
E_Hist, xedges, yedges = np.histogram2d(RRG_E[:,0], KSB_E[:,0], bins=(32,32), range = [E_Limits, E_Limits])
E_Hist.shape, xedges.shape, yedges.shape
ax.contour(xedges[0:xedges.shape[0]-1], yedges[0:yedges.shape[0]-1], E_Hist, 10, cmap = pl.cm.bone, label = 'e(1)')
ax.plot(yeqx, yeqx)
#ax.imshow(E_Hist, extent=extent, interpolation='gaussian')

ax = fE.add_subplot(2,1,2)

ax.set_xlabel('RRG Ellipticity')
ax.set_ylabel('KSB Ellipticity')

ax.set_title('e(2)')
E_Hist, xedges, yedges = np.histogram2d(RRG_E[:,1], KSB_E[:,1], bins=(32,32), range = [E_Limits, E_Limits])
E_Hist.shape, xedges.shape, yedges.shape
ax.contour(xedges[0:xedges.shape[0]-1], yedges[0:yedges.shape[0]-1], E_Hist, 10, cmap = pl.cm.bone, label = 'e(1)')
ax.plot(yeqx, yeqx)

output_name = 'Ellipticity_Comparison.pdf'
pl.savefig(output_name,format = 'pdf',bbox_inches = 'tight', pad_inches = 0.1)
print 'Produced plot output to '+output_name

pl.close()

##Produce Size Plots##
S_Input = np.genfromtxt(file_Size)
RRG_S = np.zeros((S_Input.shape[0],2))
KSB_S = np.zeros((S_Input.shape[0],2))

RRG_S[:,0] = S_Input[:,RRG_S_Col[0]]; RRG_S[:,1] = S_Input[:,RRG_S_Col[1]]
KSB_S[:,0] = S_Input[:,KSB_S_Col[0]]; KSB_S[:,1] = S_Input[:,KSB_S_Col[1]]

f = pl.figure(0, (6,2*6))
ax = f.add_subplot(2,1,1)

ax.set_xlabel('RRG Size')
ax.set_ylabel('KSB Size (Uncorrected)')

ax.set_title(r'$J_{11} + J_{22}$')
S_Hist, xedges, yedges = np.histogram2d(RRG_S[:,0], KSB_S[:,0], bins=(32,32))#, range = [S_Limits, S_Limits])
S_Hist.shape, xedges.shape, yedges.shape
ax.contour(xedges[0:xedges.shape[0]-1], yedges[0:yedges.shape[0]-1], S_Hist, 10, cmap = pl.cm.bone)

ax.set_xlim(0,10); ax.set_ylim(0,10)

ax = f.add_subplot(2,1,2)

ax.set_xlabel('RRG Size')
ax.set_ylabel('KSB Size (Uncorrected)')

ax.set_title(r'$\rm{det}(J)$')
S_Hist, xedges, yedges = np.histogram2d(RRG_S[:,1], KSB_S[:,1], bins=(32,32))#, range = [S_Limits, S_Limits])
S_Hist.shape, xedges.shape, yedges.shape
ax.contour(xedges[0:xedges.shape[0]-1], yedges[0:yedges.shape[0]-1], S_Hist, 10, cmap = pl.cm.bone)

output_name = 'Size_Comparison.pdf'
pl.savefig(output_name,format = 'pdf',bbox_inches = 'tight', pad_inches = 0.1)
print 'Produced plot output to '+output_name
            
