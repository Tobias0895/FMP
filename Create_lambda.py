# %% 
import ChiantiPy.core as ch
import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.interpolate import interp1d
import ChiantiPy.tools.data as chdata
for akey in chdata.Defaults.keys():  print(akey)
# %%
# Well try to create a lambda function fist for one wavelength
import ChiantiPy.tools.filters as chfilters

wvl = np.linspace(0.1, 100, 100) #angstrom
temperature = np.logspace(5, 9, 30)
density = 1.e+8
emeasure = np.full_like(temperature, 1.e+27)
s = ch.spectrum(temperature, density, wvl, em = emeasure, doContinuum=1, minAbund=1.e-5)
# %%

plt.plot(wvl, s.Spectrum['intensity'][15], label='T={:.2} K'.format(temperature[15]))
plt.plot(wvl, s.Spectrum['intensity'][18], label='T={:.2} K'.format(temperature[18]))
plt.xlabel(s.Spectrum['xlabel'])
plt.ylabel(s.Spectrum['ylabel'])
plt.legend()
plt.savefig('Figures/Xray spectrum.png', dpi=500)
plt.show()
# %%
# Plot lambda as a function of temperature for the whole spectrum
Intensities = np.trapz(s.Spectrum['intensity'], wvl, axis=-1)

# Interpolating this data
Linear_Interpolator = interp1d(temperature, Intensities)
quadratic_Interpolator = interp1d(temperature, Intensities, kind='quadratic')
full_temps = np.logspace(5, 9, 100)

plt.plot(full_temps, Linear_Interpolator(full_temps), label='Linear interpolated lambda')
plt.plot(full_temps, quadratic_Interpolator(full_temps), label='Quadratic interpolated lambda')
plt.scatter(temperature, Intensities, label='sampled lambda')
plt.xscale('log')
plt.xlabel('Temperature (K)')
plt.ylabel('Xray Flux (erg cm$^{-2}$ s$^{-1}$ sr$^{-1}$)')
plt.legend()
plt.savefig('Figures/Lambda_Temperature.png')
plt.show()

# %%

def create_lambda(temperatures, minwvl, maxwvl, density=1.e9, em=1.e27, continuum=1):

    wvl = np.linspace(minwvl, maxwvl, 1000) # Wavelength range for Xray
    s = ch.spectrum(temperature=temperatures, eDensity=density, wavelength=wvl, em=em, doContinuum=continuum, minAbund=1.e-4)
    Fluxes_T = np.trapz(s.Spectrum['intensity'], wvl, axis=-1) # Intergrate the spectra for each temperature to obtain the total Xray flux
    Lambda = interp1d(temperatures, Fluxes_T, kind='linear')
    return Lambda

Lambda = create_lambda(temperature, 0.1, 100)
# %%
plt.plot(full_temps, Lambda(full_temps))
plt.xscale('log')
plt.xlabel('Temperature (K)')
plt.ylabel('Xray Flux (erg cm$^{-2}$ s$^{-1}$ sr$^{-1}$)')
plt.legend()
plt.show()
    