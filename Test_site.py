# %%
from starwinds_readplt.dataset import Dataset
import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.interpolate import NearestNDInterpolator

# %%
# Import the data

ds = Dataset.from_file("/net/vdesk/data2/eikelenboom/Master/Data/data/3d__var_3_n00060000.plt")
print(ds)
# %%
# Get a stack of points and the data in a numpy format
ds_points = np.stack([ds(name) for name in ds.variables[:3]], axis=-1)
ds_data = np.stack([ds(name) for name in ds.variables], axis=-1)
#%%
# Were going to try to sum over one direction with a constant lambda

# First we create an interpolator
interpolater = NearestNDInterpolator(ds_points, ds_data)

# Create a 3d mesh grid, the different sizes is to keep track of them
x = np.linspace(-30, 30, 100)
y = np.linspace(-30, 30, 110)
z = np.linspace(-30, 30, 120)
X, Y , Z= np.meshgrid(x, y, z)
interpolated_data = interpolater(X, Y, Z)

stellar_radius = 1
mask = (X ** 2 + Y ** 2 + Z ** 2 < stellar_radius ** 2) + (X ** 2 + Y ** 2 < stellar_radius ** 2) * (Z < 0)
data_minus_star = np.where(mask==0, interpolated_data[...,ds.variables.index('Rho [g/cm^3]')], 0)
print(data_minus_star)
# %%
# integrate along the line of sight, leaving out the star and eveything behind it.
# line_of_sight_sum = np.sum(interpolated_data[...,ds.variables.index('Rho [g/cm^3]')], axis=-1)
line_of_sight_integral = np.trapz(data_minus_star, Z, axis=-1)
plt.pcolormesh(X[:, :, 0], Y[:, :, 0], line_of_sight_integral, norm='log')
# plt.colorbar()
# plt.savefig('Figures/Z_line_of_sight.pdf')
# plt.show()

# plt.pcolormesh(X[:, :, 0], Y[:, :, 0], line_of_sight_sum, norm='log')
# plt.colorbar()
# %%
def G(T):
    if T > 1e6
        return 1
    elif 
        return 0

def integrate_along_line_of_sight(direction, stellar_radius, interpolater, star_variable):

    # We first create a 3D meshgrid
    x = np.linspace(-30, 30, 100)
    y = np.linspace(-30, 30, 110)
    z = np.linspace(-30, 30, 120)
    X, Y , Z= np.meshgrid(x, y, z)

    # Interpolate the data on that grid
    interpolated_data = interpolater(X, Y, Z)
    
    # From the mesh grid create a mask that removes the star
    mask_star = X ** 2 + Y ** 2 + Z ** 2 <= stellar_radius ** 2
    # Mask the shadow of the star
    fb, car = direction.split()
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
        raise Value_error()

    density_minus_star = np.where(mask==0, interpolated_data[...,ds.variables.index('Rho [g/cm^3]')])
    temperature_minus_star = np.where(mask==0, interpolated_data[...,ds.variables.index('te [K]')])
    integrand = np.square(data)

    # Now we integrate along the line of sight
    if car.lower() == "x":
        total_flux = np.trapz(data_minus_star, X, axis=0)
    elif car.lower() == 'y':
        total_flux = np.trapz(data_minus_star, Y, axis=1)
    elif car.lower() == 'z':
        total_flux = np.trapz(data_minus_star, Z, axis=2)
    else:
        raise Value_error

    return total_flux

# %%

            
    