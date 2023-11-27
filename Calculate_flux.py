#%%
import wave
import numba
from scipy import interpolate
from starwinds_readplt.dataset import Dataset
import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.interpolate import NearestNDInterpolator, LinearNDInterpolator
from numba import jit
from tqdm import tqdm
import Rotation

__all__ = ['find_nearest', 'simple_g', 'G', 'integrate_along_line_of_sight', 'total_lum_wvl_bin',
           'create_list_of_tuples', 'create_spectra']

def import_data(loc, interpolate='nearest'):
    ds = Dataset.read()
    s = Dataset.from_file(data_loc + '/3d__var_3_n00060000.plt')

    # Get a stack of points and the data in a numpy format and create an interpolator function
    ds_points = np.stack([ds(name) for name in ds.variables[:3]], axis=-1)
    # ds_points *= solar_radius
    ds_data = np.stack([ds(name) for name in ds.variables], axis=-1)
    if interpolate.lower() == 'nearest':
        return NearestNDInterpolator(ds_points, ds_data)
    elif interpolate.lower == 'linear':
        return LinearNDInterpolator(ds_points, ds_data)
    else:
        return ds_points, ds_data



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
    low_wvl_idx = find_nearest(wvl_array, wavelength_band[0])
    high_wvl_idx = find_nearest(wvl_array, wavelength_band[1])
    G_integrated_across_band = np.trapz(G_total[...,low_wvl_idx:high_wvl_idx],
                                         wvl_array[...,low_wvl_idx:high_wvl_idx], axis=-1)
    # G_band_average = G_total[...,mean_wvl_idx]
    
    # Interpolate G for certain wavelength
    T_interpolator = interp1d(T_array, G_integrated_across_band)
    interpolated_fluxes = T_interpolator(T)
    return interpolated_fluxes


def integrate_along_line_of_sight(wvl_bin, direction, stellar_radius, interpolater, 
                                  image_radius=20, pixel_count=60, angle=0., *args, **kwargs):

    # We first create a 3D meshgrid
    image_radius *= stellar_radius
    x = np.linspace(-1*image_radius, image_radius, pixel_count)
    y = np.linspace(-1*(image_radius), image_radius, pixel_count)
    z = np.linspace(-1*image_radius, image_radius, pixel_count)
    X, Y, Z= np.meshgrid(x, y, z)
    X_rot, Y_rot, Z_rot = Rotation.rotate_grid(angle, (X, Y ,Z))
    # Interpolate the data on that grid
    interpolated_data = interpolater(X_rot, Y_rot, Z_rot)
    
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
    temperature_minus_star = np.where(mask==False, interpolated_data[...,ds.variables.index('te [K]')], 1e4)
    integrand = np.square(density_minus_star) * G(temperature_minus_star, wvl_bin, *args)
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
    # The mesh variabele is different for when we integrate from the z direction.
    if direction[1] == 'z':
        # First we integrate along one axis
        integrate_axis1 = np.trapz(LoS_flux, mesh[0], axis=-1)
        # Then integrate along the other to get a single value
        total_flux_in_band = np.trapz(integrate_axis1, mesh[0][0])
    else:
        integrate_axis1 = np.trapz(LoS_flux, mesh[1], axis=-1)
        # Then integrate along the other to get a single value
        total_flux_in_band = np.trapz(integrate_axis1, mesh[1][0])
    return total_flux_in_band


def create_list_of_tuples(lst1:list|np.ndarray, lst2:list|np.ndarray) -> list:
    result = []  # Empty list to store the tuples
    for i in range(len(lst1)):
        # Create a tuple from corresponding elements
        tuple_element = (lst1[i], lst2[i])
        result.append(tuple_element)  # Append the tuple to the list
    return result


def create_spectra(wvl_range:tuple, band_width:int|float, **kwargs) -> tuple:
    a = np.arange(wvl_range[0], wvl_range[1],band_width)
    b = np.arange(wvl_range[0] + band_width, wvl_range[1] + band_width, band_width)
    wvl_bands = create_list_of_tuples(a,b)
    center_bands = np.mean(wvl_bands, axis=-1)
    spectrum = np.array([])
    for b in wvl_bands:
        flux = total_lum_wvl_bin(b, **kwargs)
        spectrum = np.append(spectrum, flux)
    return (spectrum, center_bands)
    

if __name__ == '__main__':
    data_loc = os.environ['FMPdata']
    ds = Dataset.from_file(data_loc + '/3d__var_3_n00060000.plt')

    solar_radius = 6.934e10 # cmj0
    # Get a stack of points and the data in a numpy format and create an interpolator function
    ds_points = np.stack([ds(name) for name in ds.variables[:3]], axis=-1)
    # ds_points *= solar_radius
    ds_data = np.stack([ds(name) for name in ds.variables], axis=-1)
    mask = (ds_data[...,ds.variables.index('X [R]')] <= 10) & (ds_data[...,ds.variables.index('X [R]')] >= -10)
    reducing_mask = mask
    for i in range(ds_data.shape[1] -1):
        reducing_mask = np.column_stack((reducing_mask, mask))

    ds_reduced = np.where(reducing_mask==True, ds_data, 0)

    nearest_interpolater = NearestNDInterpolator(ds_points, ds_data)

    for angle in tqdm(range(0,360, 20)):
        band = (100, 150)
        flux, mesh = integrate_along_line_of_sight(band, direction='+x', stellar_radius=1, interpolater=nearest_interpolater,
                                                    angle=np.deg2rad(angle), image_radius=20, pixel_count=100)
        plt.pcolormesh(mesh[0], mesh[1], flux, norm='log', shading='nearest')
        plt.xlabel('Y [R]')
        plt.ylabel('R§  ')
        plt.colorbar()
        plt.axis('equal')
        # plt.savefig(f'Figures/2D_projection_{direction}_wvl_{band}')
        plt.show()