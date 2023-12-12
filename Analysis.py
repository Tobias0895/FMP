# %%
from Star_class import star_model
import numpy as np 
import matplotlib.pyplot as plt
import os
# %%
# retrieve all the stars in the data directory
names = [i for i in os.listdir(os.environ['FMPdata']) if i[0] != '.']
# make the lightcurve of the desired star
star = star_model(names[names.index('1x-Mel25-005')])
ang, lums = star.light_curve('equator')
np.save(f'light_curves/light_curve_1x-Mel25-005.npy', [ang, lums])
# %%

fig_pl = plt.figure()
ax = fig_pl.add_subplot(111)


for name in names:
    star = star_model(name)
    ax.scatter(star.params['RotationPeriodStar'] ** -2 * star.params['Radius'], np.log10(star.lum_x()), label=f'{name}')

ax.set_ylabel('Log L$_x$ (erg s$^{-1}$)')
ax.set_xlabel('P$_{rot}$ (days)')
ax.legend()
plt.show()


# %%

