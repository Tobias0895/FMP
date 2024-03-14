import ChiantiPy.core as ch
import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.interpolate import interp1d
from matplotlib.colors import LogNorm
import ChiantiPy.tools.filters as chfilters

# Well try to create a lambda function fist for one wavelength
abundance='sun_coronal_2012_schmelz_ext.abund',
wvl = np.geomspace(0.0001, 100, 3001) # angstrom
temperature = np.logspace(4, 8, 201)
abund =1e-5
density = 1.e+10
s = ch.spectrum(temperature, density, wvl, filter = (chfilters.gaussian, 0.1),
                em = 1,
                doContinuum=True,
                minAbund=abund,
                verbose=False)

L, T = np.meshgrid(wvl, temperature)
np.save('geom-wvl.npy', wvl)
np.save('temps.npy', temperature)
np.save(f'G-{abund}_geomwvl.npy', s.Spectrum['intensity'])
vmax = np.max(s.Spectrum['intensity'])
norm = LogNorm(vmax=vmax, vmin=vmax/1e6, clip=True)
plt.pcolormesh(L, T, s.Spectrum['intensity'], shading='gouraud', norm=norm, cmap='plasma')
plt.yscale('log')
plt.ylabel('log(T [K])')
plt.xlabel(s.Spectrum['xlabel'])
plt.colorbar(extend='min', label=s.Spectrum['ylabel'])
# plt.savefig('Figures/2D_G function_pretty', dpi=500)sa
plt.show()