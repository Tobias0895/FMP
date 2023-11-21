#%%
from starwinds_readplt.dataset import Dataset
import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.interpolate import NearestNDInterpolator, LinearNDInterpolator

__all__ = ['find_nearest', 'simple_g', 'G', 'integrate_along_line_of_sight', 'total_lum_wvl_bin',
           'create_list_of_tuples']

def find_nearest(array, value):
    array = np.asarray(array)
    lowest_idx = abs(array - value).argmin()
    return lowest_idx

def simple_g(T, *args):
    return np.where((T >=1e7) * (T <= 1e8), 1, 0)

def G(T, wavelength_band:tuple):
    from scipy.interpolate import interp1d
    # Load in the created G function and the corresponding wvl,T grid
    G_total = np.load('G-0.0001.npy')
    wvl_array = np.load('L.npy')[0]
    T_array = np.load('T.npy')[...,0]
    
    # Find the indicies of the closest wavelengths on the grid
    mean_wvl = np.mean(wavelength_band)
    mean_wvl_idx = find_nearest(wvl_array, mean_wvl)
    G_band_average = G_total[...,mean_wvl_idx]
    
    # Interpolate G for certain wavelength
    T_interpolator = interp1d(T_array, G_band_average)
    interpolated_fluxes = T_interpolator(T)
    return interpolated_fluxes


def integrate_along_line_of_sight(wvl_bin, direction, stellar_radius, interpolater, image_radius=20, pixel_count=60, *args, **kwargs):

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
    density_minus_star = np.where(mask==False, interpolated_data[...,ds.variables.index('Rho [g/cm^3]')], 0)
    temperature_minus_star = np.where(mask==False, interpolated_data[...,ds.variables.index('te [K]')], 0)
    integrand = density_minus_star #* G(temperature_minus_star, wvl_bin, *args)
    # integrand = np.square(density_minus_star) * simple_g(temperature_minus_star, wvl_bin, *args)

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
    
def total_lum_wvl_bin(wvl_bin:tuple, direction:str, stellar_radius, interpolater, *args, **kwargs) -> float:
    LoS_flux, mesh = integrate_along_line_of_sight(wvl_bin, direction, stellar_radius, interpolater=interpolater, *args, **kwargs)
    # First we integrate along one axis
    integrate_axis1 = np.trapz(LoS_flux, mesh[0], axis=1)
    # Then integrate along the other to get a single value
    total_flux_in_band = np.trapz(integrate_axis1, mesh[0][0])
    return total_flux_in_band


def create_list_of_tuples(lst1, lst2):
    result = []  # Empty list to store the tuples
    for i in range(len(lst1)):
        # Create a tuple from corresponding elements
        tuple_element = (lst1[i], lst2[i])
        result.append(tuple_element)  # Append the tuple to the list
    return result


data_loc = os.environ['FMPdata']
ds = Dataset.from_file(data_loc + '/3d__var_3_n00060000.plt')

# Get a stack of points and the data in a numpy format and create an interpolator function
ds_points = np.stack([ds(name) for name in ds.variables[:3]], axis=-1)
ds_data = np.stack([ds(name) for name in ds.variables], axis=-1)
mask = (ds_data[...,ds.variables.index('X [R]')] <= 10) & (ds_data[...,ds.variables.index('X [R]')] >= -10)
reducing_mask = mask
for i in range(ds_data.shape[1] -1):
    reducing_mask = np.column_stack((reducing_mask, mask))
print(reducing_mask.shape)

ds_reduced = np.where(reducing_mask==True, ds_data, 0)

print(ds_reduced[...,:3])
nearest_interpolater = NearestNDInterpolator(ds_reduced[...,:3], ds_reduced)
# linear_interpolater = LinearNDInterpolator(ds_points[reducing_mask], ds_reduced)


# We create wavelength bins 
band_width=3
a = np.arange(0,100,band_width)
b = np.arange(band_width,100 + band_width,band_width)
wavelength_band = create_list_of_tuples(a,b)
center_bands = np.mean(wavelength_band, axis=-1)
# %%
band_fluxes = []
for band in wavelength_band:
    total_band_flux = total_lum_wvl_bin(band, direction='+z', stellar_radius=1, interpolater=nearest_interpolater_interpolater, pixel_count=100)
    band_fluxes.append(total_band_flux)

total_flux = np.trapz(band_fluxes, center_bands)
print(total_flux)
bars = plt.bar(center_bands, band_fluxes, width=band_width, log=True, edgecolor='k', align='center')
plt.ylabel('erg $cm^-2$ s$^-1$ ')
plt.xlabel('Wavelength ($\AA$)')
plt.savefig("Figures/Intensity_+z.png", dpi=500)
    # %%
flux, mesh = integrate_along_line_of_sight((16, 18), direction='+z', stellar_radius=1, interpolater=nearest_interpolater, image_radius=60, pixel_count=100)
print(np.max(flux))
plt.pcolormesh(mesh[0], mesh[1], flux, norm='log')
plt.xlabel('X [R]')
plt.ylabel('Y [R]')
plt.colorbar()
plt.axis('equal')
plt.show()