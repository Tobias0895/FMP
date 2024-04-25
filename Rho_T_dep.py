from Xray_winds.src.Xray_winds import load_data
from Xray_winds.src.Xray_winds import Star_class
from Xray_winds.src.Xray_winds import Grid_Operations as go
import numpy as np
from matplotlib import pyplot as plt
import os
import pandas as pd
import astropy.constants as c
from pyvo import registry
from tqdm import tqdm


mnras_single_width_in = 3.3209

model_names = [i for i in os.listdir(os.environ['FMPdata']) if i[0] != '.']
star_data = pd.read_csv('/Users/tobias/Data/Extended_stellar_params.csv', dtype=None, delimiter=';', encoding='utf-8')

fig = plt.figure(figsize=(6,6))
figTR = plt.figure(figsize=(6,6))
axTR = figTR.add_subplot(111)
ax = fig.add_subplot(121)
ax2 = fig.add_subplot(122)


for name in tqdm(model_names):
    if "sun" in name.lower(): # Sun data is not found in star_data so we make an exception for those 
        association = 'Sun'
        rossby = 2
    else:
        star_idx= star_data[star_data['Star'].str.replace(' ', '') == name.replace('1x-', '')].index
        association = star_data['Associations'][star_idx].values[0]

    star = Star_class.star_model(name, interpolation='nearest')
    ds, *_ = star.raw_data()
    data = np.stack([ds(name) for name in ds.variables], axis=-1)
    points = np.stack([ds(name) for name in ds.variables[:3]], axis=-1)

    temperature = data[...,star.var_list.index('te [K]')]
    density = data[...,star.var_list.index('Rho [g/cm^3]')]
    log_d = np.log10(density)
    log_T = np.log10(temperature)
    r = np.sqrt(points[...,0] ** 2 + points[...,1] ** 2 + points[...,2] ** 2)

    unique_r = np.unique(r)
    meanr_T = []
    meanr_d = []
    for rn in unique_r:
        radial_mask = r == rn
        meanr_T.append(np.nanmean(temperature[radial_mask]))
        meanr_d.append(np.nanmean(log_d[radial_mask]))

    maxT_radius = unique_r[np.argmax(meanr_T)]
    (Xs, Ys, Zs) = go.fibonacci_sphere(maxT_radius * star.params['RadiusStar'], 50)
    interpolated_T = star.interpolator(Xs, Ys, Zs)[...,star.var_list.index('te [K]')]
    interpolated_d = star.interpolator(Xs, Ys, Zs)[...,star.var_list.index('Rho [g/cm^3]')]
    kpr= np.log10(1.86e-3 * star.params['RotationPeriodStar'] ** -2 * (star.params['RadiusStar'] / c.R_sun.to('cm').value) **-4)
    axTR.scatter(kpr, (np.mean(interpolated_d) ** 2) * np.mean(interpolated_T), s=100, marker='*', label=association)


handles, labels = axTR.get_legend_handles_labels()
by_label = dict(zip(labels, handles))
axTR.legend(by_label.values(), by_label.keys())
axTR.invert_xaxis()
axTR.set_xlabel('k P$^{-2}$ R$^{-4}$')
axTR.set_ylabel('T$_{e}$ [K]')
ax.set_xscale('log') ; ax2.set_xscale('log')
ax.set_xlabel('Radius'); ax2.set_xlabel('Radius')
ax.set_ylabel('Log10(N$_{h}$ [g/cm$^3$])'); ax2.set_ylabel('T$_e$ [K]')
ax.set_title('Density') ; ax2.set_title('Temperature')
ax2.legend()
# fig.savefig('Figures/RhoandT_plot.png', dpi=100, bbox_inches='tight')
# figTR.savefig('Figures/MeanTempVRadius.pdf', dpi=100, bbox_inches='tight')
plt.show() 