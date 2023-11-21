# %%
import ChiantiPy.core as ch
import ChiantiPy
import numpy as np
import matplotlib.pyplot as plt

# %%
fe14 = ch.ion('fe_14',temperature=2e6, eDensity=1e9)

# wvl = 1 + 0.01 * np.arange(100)
wvl = 0.5 + 0.02 * np.arange(8001)
fe14.spectrum(wvl)

c = ch.continuum('fe_14', temperature = 0.5e6, em=1e+27)
c.freeFree(wvl)
c.freeBound(wvl)
plt.plot(wvl, c.FreeFree['intensity'], label='ff')
# plt.plot(wvl, fe14.Spectrum["intensity"])
plt.plot(wvl, c.FreeBound['intensity'], label='fb')
plt.plot(wvl, c.FreeFree['intensity'] + c.FreeBound['intensity'])
plt.xlabel(c.FreeFree['xlabel'])
plt.ylabel(c.FreeFree['ylabel'])
plt.legend()
plt.show()

# %%
# Try the same with hydrogen

h1 = ch.ion('h_1', temperature=283400, eDensity=1.0877e-19 * 6.02e22)
print(h1.Spectroscopic)
wvl2 = np.linspace(0.1, 100, 1000)
h1.spectrum(wvl2)

c1 = ch.continuum('h_1', temperature=283400)
c.freeFree(wvl2)
c.freeBound(wvl2)
plt.plot(wvl2, c.FreeFree['intensity'], label='FreeFree')
plt.plot(wvl2, c.FreeBound['intensity'], label='FreeBound')
plt.plot(wvl2, c.FreeBound['intensity']+ c.FreeFree['intensity'], label='FreeFree + FreeBound')
plt.plot(wvl2, h1.Spectrum['intensity'], l  abel='Lines')
plt.xlabel(h1.Spectrum['xlabel'])
plt.ylabel(h1.Spectrum['ylabel'])
plt.legend()
plt.show()
# %%
# Multi-Bunch Ion class

Xwvl = np.linspace(0.1, 100, 1000)
Abundicies = [1e-4]
temperature = 1e6 #np.linspace(1e5, 1e7, 100)
density = 1.e9
for ab in Abundicies:
    s = ch.spectrum(temperature, density, Xwvl, doContinuum=1, minAbund=ab)
    plt.plot(Xwvl, s.Spectrum['intensity'], label=f'min {ab}')
plt.legend()
plt.xscale('log')
plt.show()
