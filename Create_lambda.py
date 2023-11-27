# %% 
import ChiantiPy.core as ch
import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.interpolate import interp1d

# %%
# Well try to create a lambda function fist for one wavelength
import ChiantiPy.tools.filters as chfilters
abundance='sun_coronal_2012_schmelz_ext.abund',
wvl = np.linspace(0.1, 180, 2001) # angstrom
temperature = np.logspace(4, 8, 201)
abund =2e-5
density = 1.e+10
s = ch.spectrum(temperature, density, wvl, 
                filter = (chfilters.gaussian, 0.1),
                em = None,
                doContinuum=True,
                minAbund=abund,
                abundance=abundance,
                verbose=False)
# %%
from matplotlib.colors import LogNorm

L, T = np.meshgrid(wvl, temperature)
# np.save('L.npy', L)
# np.save('T.npy', T)
# np.save(f'G-{abund}.npy', s.Spectrum['intensity'])
vmax = np.max(s.Spectrum['intensity'])
norm = LogNorm(vmax=vmax, vmin=vmax/1e6, clip=True)
plt.pcolormesh(L, T, s.Spectrum['intensity'], shading='gouraud', norm=norm, cmap='plasma')
plt.yscale('log')
plt.ylabel('log(T [K])')
plt.xlabel(s.Spectrum['xlabel'])
plt.colorbar(extend='min', label=s.Spectrum['ylabel'])
# plt.savefig('Figures/2D_G function_pretty', dpi=500)sa
plt.show()
# %%

    