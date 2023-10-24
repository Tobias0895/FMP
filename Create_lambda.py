# %% 
import ChiantiPy.core as ch
import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.interpolate import LinearNDInterpolator
import ChiantiPy.tools.data as chdata
for akey in chdata.Defaults.keys():  print(akey)
# %%
# Well try to create a lambda function fist for one wavelength
import ChiantiPy.tools.filters as chfilters

wvl = np.linspace(0.1, 100, 100) #angstrom
temperature = np.logspace(2e5, 2e9, 30)
density = 1.e+9
emeasure = np.full_like(temperature, 1.e+27)
s = ch.spectrum(temperature, density, wvl, filter = (chfilters.gaussian,.2), em = emeasure, doContinuum=0, minAbund=1.e-5)
print(s)
