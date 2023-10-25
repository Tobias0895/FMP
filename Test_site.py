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
x = np.linspace(-30, 30, 60)
y = np.linspace(-30, 30, 60)
z = np.linspace(-30, 30, 60)
X, Y , Z= np.meshgrid(x, y, z)
interpolated_data = interpolater(X, Y, Z)

stellar_radius = 8
mask = (X ** 2 + Y ** 2 + Z ** 2 < stellar_radius ** 2) + (Z ** 2 + Y ** 2 < stellar_radius ** 2) * (X > 0)
data_minus_star = np.where(mask==0, interpolated_data[...,ds.variables.index('Rho [g/cm^3]')], 0)
print(data_minus_star)
# %%
ax = plt.figure().add_subplot(projection='3d')
ax.voxels(mask)
ax.axis('equal')
plt.savefig('Figures/Masked_cubes.png', dpi=500)