# %%
from starwinds_readplt.dataset import Dataset
import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.interpolate import NearestNDInterpolator
# import numba

# %%
# @jit(nonpython=True)
# Import the data
data_loc = os.environ['FMPdata']
ds = Dataset.from_file(data_loc + '/3d__var_3_n00060000.plt')
print(ds)
# %%
# Get a stack of points and the data in a numpy format
ds_points = np.stack([ds(name) for name in ds.variables[:3]], axis=-1)
ds_data = np.stack([ds(name) for name in ds.variables], axis=-1)
interpolater = NearestNDInterpolator(ds_points, ds_data)
# %%
def find_nearest(array, value):
    array = np.asarray(array)
    lowest_idx = abs(array - value).argmin()
    return lowest_idx


def G(T, wavelength_band:tuple):
    G_total = np.load('G-5e-05.npy')
    wvl_array = np.load('L.npy')[0]
    T_array = np.load('T.npy')[0]
    min_wvl_idx = find_nearest(wvl_array, wavelength_band[0])
    max_wvl_idx = find_nearest(wvl_array, wavelength_band[1])
    G_band = G_total[...,min_wvl_idx:max_wvl_idx]




    mask = T > 1e6
    G_temps = np.where(mask==0, T, 0) # where the mask is False keep T, else 0
    return G_temps

def integrate_along_line_of_sight(direction, stellar_radius, interpolater, image_radius=40, pixel_count=60, *args, **kwargs):

    # We first create a 3D meshgrid
    x = np.linspace(-1*image_radius, image_radius, pixel_count)
    y = np.linspace(-1*image_radius, image_radius, pixel_count)
    z = np.linspace(-1*image_radius, image_radius, pixel_count)
    X, Y, Z= np.meshgrid(x, y, z)

    # Interpolate the data on that grid
    interpolated_data = interpolater(X, Y, Z)
    
    # From the mesh grid create a mask that removes the star
    mask_star = X ** 2 + Y ** 2 + Z ** 2 <= stellar_radius ** 2
    # Mask the shadow of the star
    fb = list(direction)[0]
    car = list(direction)[1]
    if car.lower() == 'x':
        mask_shadow = Y ** 2 + Z ** 2 <= stellar_radius ** 2
        if fb == "+":
           mask_shadow *= (X < 0)
        elif fb == '-':
            mask_shadow *= (X > 0)
    elif car.lower() == 'y':
        mask_shadow = X ** 2 + Z ** 2 <= stellar_radius ** 2
        if fb == "+":
            mask_shadow *= (Y < 0)
        elif fb == '-':
            mask_shadow *= (Y > 0)
    elif car.lower() == 'z':
        mask_shadow = X ** 2 + Y ** 2 <= stellar_radius ** 2
        if fb == "+":
            mask_shadow *= (Z < 0)
        elif fb == '-':
            mask_shadow *= (Z > 0)   
    else:
        raise ValueError
    mask = mask_star + mask_shadow
    density_minus_star = np.where(mask==0, interpolated_data[...,ds.variables.index('Rho [g/cm^3]')], 0)
    temperature_minus_star = np.where(mask==0, interpolated_data[...,ds.variables.index('te [K]')], 0)
    integrand = np.square(density_minus_star) * G(temperature_minus_star, *args)

    # Now we integrate along the line of sight
    if car.lower() == "x":
        total_flux = np.trapz(integrand, X, axis=1)
        return total_flux, (Y[:, 0, :], Z[:, 0, :])
    elif car.lower() == 'y':
        total_flux = np.trapz(integrand, Y, axis=0)
        return total_flux, (X[0,...], Z[0,...])
    elif car.lower() == 'z':
        total_flux = np.trapz(integrand, Z, axis=-1)
        return total_flux, (X[...,0], Y[...,0])
    else:
        raise ValueError
    
def total_lum_wvl_bin(wvl_bin:tuple, **kwargs) -> float:
    LoS_flux, mesh = integrate_along_line_of_sight(wvl_bin, **kwargs)
    Los_array = np.ravel(LoS_flux)
    total_lux_in_band = np.trapz(Los_array, dx=np.diff(mesh[0])[0,0])
    return total_lux_in_band



# %%
tobi, mesh = integrate_along_line_of_sight('+z', 1, interpolater, pixel_count=100)
# print(mesh.shape)
plt.pcolormesh(mesh[0], mesh[1], tobi, norm='log', shading='gouraud')
plt.axis('equal')
plt.colorbar()
plt.savefig('Figures/intergrated_line_of_sight_x.png', dpi=500)
plt.show()

total_lum_wvl_bin((0,2), direction='+z', stellar_radius=1, interpolater=interpolater, pixel_count=50)
# %%
wavelength_band = (6,8)
G_total = np.load('G-5e-05.npy')
print(G_total.shape)
wvl_array = np.load('L.npy')[0]
T_array = np.load('T.npy')[0]
a = find_nearest(wvl_array, wavelength_band[0])
b = find_nearest(wvl_array, wavelength_band[1])
G_total[...,a:b].shape

