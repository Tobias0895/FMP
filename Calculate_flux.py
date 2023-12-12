
from venv import create
import numpy as np
import matplotlib.pyplot as plt
from sympy import fraction
from tqdm import tqdm
import Rotation

__all__ = ['find_nearest', 'simple_g', 'G', 'projection_2d', 'total_lum_wvl_bin',
           'create_list_of_tuples', 'create_spectra']

def find_nearest(array, value):
    array = np.asarray(array)
    lowest_idx = abs(array - value).argmin()
    return lowest_idx

def simple_g(T, *args):
    return np.where((T >=1e6) * (T <= 1e7), 1, 0)

def G(T, wavelength_band:tuple, use_simple_g=False):
    from scipy.interpolate import interp1d
    if use_simple_g == True:
        return simple_g(T)
    else:
        # Load in the created G function and the corresponding wvl,T grid
        G_total = np.load('G_(T,L)/G-0.0001.npy')
        wvl_array = np.load('G_(T,L)/L.npy')[0]
        T_array = np.load('G_(T,L)/T.npy')[...,0]
        
        # Find the indicies of the closest wavelengths on the grid
        low_wvl_idx = find_nearest(wvl_array, wavelength_band[0])
        high_wvl_idx = find_nearest(wvl_array, wavelength_band[1])
        G_integrated_across_band = np.trapz(G_total[...,low_wvl_idx:high_wvl_idx],
                                            wvl_array[...,low_wvl_idx:high_wvl_idx], axis=-1)
        
        # Interpolate G for certain wavelength
        T_interpolator = interp1d(T_array, G_integrated_across_band)
        interpolated_fluxes = T_interpolator(T)
        return interpolated_fluxes


def projection_2d(wvl_bin:tuple, stellar_radius:float, interpolator, var_list:list,
                                  image_radius=20, pixel_count=60, angle=(0.,0.), direction='+x', *args, **kwargs) -> tuple:
    # We first create a 3D meshgrid
    image_radius *= stellar_radius
    x = np.linspace(-1*image_radius, image_radius, pixel_count)
    y = np.linspace(-1*(image_radius), image_radius, pixel_count)
    z = np.linspace(-1*image_radius, image_radius, pixel_count)
    X, Y, Z= np.meshgrid(x, y, z)
    X_rot, Y_rot, Z_rot = Rotation.rotate_grid(angle[0], angle[1], (X, Y ,Z))
    # Interpolate the data on that grid
    interpolated_data = interpolator(X_rot, Y_rot, Z_rot)

    
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
    # density_minus_star = np.where(mask==False, interpolated_data[...,var_list.index('Rho [g/cm^3]')], 0)
    # temperature_minus_star = np.where(mask==False, interpolated_data[...,var_list.index('te [K]')], np.nan)
    integrand = np.square(interpolated_data[...,var_list.index('Rho [g/cm^3]')] / 1.67e-24) * G(interpolated_data[...,var_list.index('te [K]')], wvl_bin, *args, **kwargs)
    masked_integrand = np.where(mask==False, integrand, 0)

    
    # We want to take into account that light that goes towards the star can't be accounted for
    Area_of_star = np.square(stellar_radius) * np.pi 
    # Creating an array in the shape of the data to calculate the solid angle of the star at each point
    solid_angle_array = Area_of_star / (X**2 + Y**2 + Z**2)
    # The fraction of the sky at a point is then solid_angle / 4pi
    fraction_usable_light =  1 - (solid_angle_array / (4*np.pi))
    masked_integrand *= fraction_usable_light

    # Now we integrate along the line of sight
    if car.lower() == "x":
        total_flux = np.trapz(masked_integrand, X, axis=1)
        return total_flux, (Y[:, 0, :], Z[:, 0, :])
    elif car.lower() == 'y':
        total_flux = np.trapz(masked_integrand, Y, axis=0)
        return total_flux, (X[0,...], Z[0,...])
    elif car.lower() == 'z':
        total_flux = np.trapz(masked_integrand, Z, axis=-1)
        return total_flux, (X[...,0], Y[...,0])
    else:
        raise ValueError
    
def total_lum_wvl_bin(wvl_bin:tuple, stellar_radius, interpolator, *args, **kwargs) -> float:
    LoS_flux, mesh = projection_2d(wvl_bin, stellar_radius, interpolator=interpolator, *args, **kwargs)
    # The mesh variabele is different for when we integrate from the z direction.
    # if direction[1] == 'z':
    #     # First we integrate along one axis
    #     integrate_axis1 = np.trapz(LoS_flux, mesh[0], axis=-1)
    #     # Then integrate along the other to get a single value
    #     total_flux_in_band = np.trapz(integrate_axis1, mesh[0][0])
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


def create_spectra(wvl_range:tuple, band_width:int|float, disable_tqdm=False, save_spectra=False, **kwargs) -> tuple:
    """_summary_

    Args:
        wvl_range (tuple): _description_
        band_width (int | float): _description_
        disable_tqdm (bool, optional): _description_. Defaults to False.
        save_spectra (str, optional): If True saves the spectra. Defaults to False.

    Returns:
        tuple: _description_
    """    ''''''
    a = np.arange(wvl_range[0], wvl_range[1],band_width)
    b = np.arange(wvl_range[0] + band_width, wvl_range[1] + band_width, band_width)
    wvl_bands = create_list_of_tuples(a,b)
    center_bands = np.mean(wvl_bands, axis=-1)
    spectrum = np.array([])
    for b in tqdm(wvl_bands, disable=disable_tqdm):
        flux = total_lum_wvl_bin(b, **kwargs)
        spectrum = np.append(spectrum, flux)
    if save_spectra:
        np.save(save_spectra, [center_bands, spectrum]) # type: ignore
        return (center_bands, spectrum)
    else:
        return (center_bands, spectrum)