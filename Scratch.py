#%%
import numpy as np
import ChiantiPy.core as ch
import matplotlib.pyplot as plt
from ChiantiPy.tools import filters
from matplotlib.colors import LogNorm

radio_freq = np.linspace(3, 300e9) # Hz
radio_wvl = 2.998e18 / radio_freq #AA
temps = 1e7 #np.logspace(4, 8, 201)
density = 1.e10
s = ch.spectrum(temps, density, radio_wvl,
                filter = (filters.gaussian, 0.1),
                doContinuum=True,
                verbose=False,
                minAbund=1.e-6)
# %%
L, T = np.meshgrid(radio_wvl, temps)
vmax = np.max(s.FreeFree['intensity'])
norm = LogNorm(vmax=vmax, vmin=vmax/1e6, clip=True)
fig = plt.figure()
ax = fig.add_subplot(121)
# mesh = ax.pcolormesh(L ,T, s.Spectrum['intensity'], norm=norm)
ax.plot(radio_wvl, s.Spectrum['intensity'])
plt.show()