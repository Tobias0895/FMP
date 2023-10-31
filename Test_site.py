# %%
from starwinds_readplt.dataset import Dataset
import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.interpolate import NearestNDInterpolator

# %%
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
def G(T):
    mask = T <= 1e6
    G_temps = np.where(mask==0, T, 0)
    return G_temps

def integrate_along_line_of_sight(direction, stellar_radius, interpolater):

    # We first create a 3D meshgrid
    x = np.linspace(-20, 20, 60)
    y = np.linspace(-20, 20, 60)
    z = np.linspace(-20, 20, 60)
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
    integrand = np.square(density_minus_star) * G(temperature_minus_star)

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

# %%
tobi, mesh = integrate_along_line_of_sight('+y', 1, interpolater)
plt.pcolormesh(mesh[0], mesh[1], tobi, norm='log')
plt.axis('equal')
plt.colorbar()
plt.show()