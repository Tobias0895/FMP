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
y = np.linspace(-30, 30, 61)
z = np.linspace(-30, 30, 62)
X, Y, Z = np.meshgrid(x, y, z)
interpolated_data = interpolater(X, Y, Z)
print(interpolated_data[...,ds.variables.index("Rho [g/cm^3]")].shape)

z_sum = np.sum(interpolated_data[...,ds.variables.index("Rho [g/cm^3]")], axis=-1)
print(z_sum.shape)

# %%
print(X[:, :, 0], X[:, :, 0].shape)
plt.pcolormesh(X[:, :, 0], Y[:, :, 0], z_sum, norm='log')
plt.colorbar()
plt.show()

plt.pcolormesh(X[:, :, 0], Y[:, :, 0], interpolated_data[...,ds.variables.index("Rho [g/cm^3]")][:, :, 31], norm='log')
plt.colorbar()
plt.show()
# %%
# The star is just the central cube, so we need to set everything behind that at 0
total = np.array([])
densities = interpolated_data[...,ds.variables.index('Rho [g/cm^3]')]
for i in range(densities.shape[-1]):
    sl = densities[...,i]
    