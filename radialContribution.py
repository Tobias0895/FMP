import Xray_winds.src.Xray_winds.load_data as ld
from Xray_winds.src.Xray_winds.Calculate_flux import G
import Xray_winds.src.Xray_winds.Grid_Operations as go
import scipy.interpolate as interpolate
import matplotlib.pyplot as plt
import os
import numpy as np



names = [i for i in os.listdir(os.environ['FMPdata']) if i[0] != '.']
contributions_star = []
radiis = np.linspace(1, 250, 50)
for name in names:
    interpolator, var_list, params, _, ds_points, _ = ld.import_data(name, interpolate='nearest', full_output=True)
    radii = radiis * params['RadiusStar']
    dr = np.diff(radii)[0]

    (X, Y, Z), x = go.create_grid(250 * params['RadiusStar'], 200, type='linear')
    mask_star = X ** 2 + Y ** 2 + Z ** 2 <= params['RadiusStar'] ** 2

    interpolated_data = interpolator(X, Y ,Z)

    integrand = np.square(interpolated_data[...,var_list.index('Rho [g/cm^3]')] / 1.67e-24) * G(interpolated_data[...,var_list.index('te [K]')], (0.1, 180))
    masked_integrand = np.where(mask_star==False, integrand, 0)
    contributions = []
    contributions_star.append(contributions)

    for r in radii:
        u,v,w = go.fibonacci_sphere(r, 500)

        tobi = interpolate.interpn((x,x,x), masked_integrand, (u,v,w))
        contributions.append(np.mean(tobi) * 4 * np.pi * np.square(r) * dr) # the 4pir^2 dr to account for the volume added


    
    total_contribution = np.trapz(contributions, radii)

con = np.log10(contributions_star)

mean_con = np.mean(con, axis=0)
std_con = np.std(con, axis=0)
fig, ax = plt.subplots(1,1, figsize=(8,8))
ax.plot(radiis, mean_con, label='Mean Contribution', c='k')
ax.fill_between(radiis, mean_con - std_con, mean_con + std_con, alpha=0.2, color='k', label='1 Std contribution')
ax.legend()
ax.set_xlabel('Radius (r$_{star}$)')
ax.set_ylabel('Log Lx (erg s$^{}$)')

plt.savefig('mean_radial_contributions.pdf', dpi=100, bbox_inches='tight')
plt.show()