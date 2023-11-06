# %% 
import ChiantiPy.core as ch
import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.interpolate import interp1d
# %%
# Well try to create a lambda function fist for one wavelength
import ChiantiPy.tools.filters as chfilters

wvl = np.geomspace(0.1, 50, 5000) #angstrom
temperature = np.logspace(4, 8, 100)

density = 1.e+9
emeasure = np.full_like(temperature, 1e27)
s = ch.spectrum(temperature, density, wvl, filter = (chfilters.gaussian, 0.5), em = emeasure, doContinuum=1, minAbund=1.e-5, verbose=0)
# %%
import matplotlib as mpl
L, T = np.meshgrid(wvl, temperature[60:])

plt.pcolormesh(L, T, s.Spectrum['intensity'][60:], shading='gouraud', norm='log')
plt.yscale('log')
plt.ylabel('log(T [K])')
plt.xlabel(s.Spectrum['xlabel'])
plt.savefig('Figures/2D_G function', dpi=500)
plt.colorbar()
plt.show()
# %%

    