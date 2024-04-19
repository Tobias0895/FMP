from scipy import interpolate
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

CATALOGUE = "J/ApJ/743/48"
catalogue_ivoid = f"ivo://CDS.VizieR/{CATALOGUE}"
voresource = registry.search(ivoid=catalogue_ivoid)[0]
tables = voresource.get_tables()
table_name = "J/ApJ/743/48/catalog"

tap_service = voresource.get_service("tap")
tap_records = voresource.get_service("tap").run_sync(
    f'select  * from "{table_name}"',
)
 
dfw = tap_records.table.to_pandas()
print(dfw.keys())


model_names = [i for i in os.listdir(os.environ['FMPdata']) if i[0] != '.']
star_data = pd.read_csv('/Users/tobias/Data/Extended_stellar_params.csv', dtype=None, delimiter=';', encoding='utf-8')
(X, Y, Z), _ = go.create_grid(250, 50, 'linear')
 
fig = plt.figure(figsize=(16,8))
figTR = plt.figure()
axTR = figTR.add_subplot(111)
ax = fig.add_subplot(121)
ax2 = fig.add_subplot(122)


for name in tqdm(model_names):
    
    star = Star_class.star_model(name, interpolation='nearest', verbose=False)
    ds, *_ = star.raw_data()
    data = np.stack([ds(name) for name in ds.variables], axis=-1)
    points = np.stack([ds(name) for name in ds.variables[:3]], axis=-1)

    
    temperature = data[...,star.var_list.index('te [K]')]
    density = data[...,star.var_list.index('Rho [g/cm^3]')]
    log_d = np.log10(density)
    log_T = np.log10(temperature)
    r = np.sqrt(points[...,0] ** 2 + points[...,1] ** 2 + points[...,2] ** 2)
    unique_r = np.unique(r)
    unique_T = []
    for rn in unique_r:
        unique_T.append(np.nanmean(unique_T))

 
    kpr_value = np.log10(1.86e-3 * star.params['RotationPeriodStar'] ** -2 * (star.params['RadiusStar'] / c.R_sun.to('cm').value) **-4)
    # ax.errorbar(kpr_value, np.mean(log_d), yerr=0.1 * np.std(log_d), c='k', fmt='.', markersize=5, label='raw')
    # ax2.errorbar(kpr_value, np.mean(log_T), yerr=0.1 * np.std(log_T), c='r', fmt='.', markersize=5, label='raw')
    axTR.scatter(unique_r, unique_T, s=2)

    # interp_T = np.log10(interp_data[...,star.var_list.index('te [K]')])
    # interp_d = np.log10(interp_data[...,star.var_list.index('Rho [g/cm^3]')])
    # ax.errorbar(kpr_value, np.mean(interp_d), yerr=0.1 * np.std(interp_d), label='Interpolated')
    # ax2.errorbar(kpr_value, np.mean(interp_T), yerr=0.1 * np.std(interp_T), label='Interpolated')
axTR.set_xscale('log')
ax.invert_xaxis(); ax2.invert_xaxis()
ax.set_xlabel('Kpr'); ax2.set_xlabel('KPR')
ax.set_ylabel('Log10(N$_{h}$ [g/cm$^3$])'); ax2.set_ylabel('Log10(T$_e$ [K])')
# ax2.scatter(np.log10(1.86e-3 * dfw['R*']**-4 *dfw['Prot'] **-2), dfw['Lx/bol'], c='gray', s=2, label='Reiners+2014')
# ax.scatter(np.log10(1.86e-3 * dfw['R*']**-4 *dfw['Prot'] **-2), dfw['Lx/bol'], c='gray', s=2, label='Reiners+2014')
# ax.legend(); ax2.legend()
ax.set_title('Density') ; ax2.set_title('Temperature')
fig.savefig('Figures/RhoandT_plot.png', dpi=100, bbox_inches='tight')
figTR.savefig('Figures/radiusVtemperature.png')
plt.show()