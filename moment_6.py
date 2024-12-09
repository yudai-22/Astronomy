from astropy.io import fits
import numpy as np
from matplotlib import pyplot as plt
import aplpy
import pylab
import matplotlib

fname = ".fits"
fname_new = "mom6.fits"

hdu = fits.open(fname)[0]  #[0] => primary HDU(header and data)

#make mom6
#specify the range of channels to use
start_ch = 0
end_ch = 81

ndata = hdu.data.copy()
squared_ndata = np.square(ndata[start_ch:end_ch])
mean_ndata = np.mean(squared_ndata,axis=0)
rms_ndata = np.sqrt(mean_ndata)

#make mom6 HDU
new_data = rms_ndata
new_header = hdu.header.copy()

new_header.pop('C*3')
new_header.pop('NAXIS3')
new_header.pop('PC*')
new_header.update(NAXIS=2)
nhdu = fits.PrimaryHDU(new_data,new_header)
nhdu.writeto( fname_new , overwrite=True) 

#make mom6 pdf
pylab.rcParams['font.family'] = 'serif'
pylab.rcParams['lines.linewidth'] = 0.5
matplotlib.rcParams["font.family"] = "serif"
matplotlib.rcParams["font.size"] = 15
pylab.rcParams["xtick.direction"] = "in"
pylab.rcParams["ytick.direction"] = "in"

fig=aplpy.FITSFigure(fname_new)
fig.show_colorscale(vmin=0.04,vmax=0.12,cmap='nipy_spectral',stretch='linear')
fig.add_colorbar()
fig.add_beam(color='k')
fig.colorbar.set_axis_label_text("[ K ]")
fig.colorbar.set_axis_label_font(size=15)
fig.add_scalebar(length=(2000.0/3600)/43,label='2$\,$kpc',linewidth=1.,color='black',corner="bottom left")
fig.recenter(143.042,21.504,width=0.0715103,height=0.0943856)  #degrees
fig.tick_labels.set_font(size=12) 
fig.ticks.set_color("black")
fig.add_label(0.15,0.95,"",relative=True,color='black',family="serif",size=20)

fig.savefig('_mon6_python.pdf',dpi=200)


