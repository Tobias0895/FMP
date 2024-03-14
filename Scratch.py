#%%
from networkx import project
import numpy as np
import ChiantiPy.core as ch
import matplotlib.pyplot as plt
from ChiantiPy.tools import filters
from matplotlib.colors import LogNorm
import random

# %%
import os

from pyparsing import delimitedList
data_loc = os.environ['FMPdata']
model_name = '1x-PWAnd'
star_file = open(os.path.join(data_loc, model_name + '/STAR.in'), 'r')
param_lines = star_file.readlines()
star_params = {}
for n, line in enumerate(param_lines[-3:]):
    split = line.split('  ')
    value = split[0]
    name = split[-1].split(' ')
    if name[0].isalpha() == False:
        name = split[-1].split(' ')[1]
    else:
        name = name[0]
    star_params[str(name)] = value
print(star_params['RotationPeriodStar'])


# %%
import argparse

def fibonacci_sphere(num_points: int, r: float):
    ga = (3 - np.sqrt(5)) * np.pi # golden angle                                                                             

    # Create a list of golden angle increments along tha range of number of points                                           
    theta = ga * np.arange(num_points)

    # Z is a split into a range of -1 to 1 in order to create a unit circle                                                  
    z = np.linspace(1/num_points-1, 1-1/num_points, num_points)

    # a list of the radii at each height step of the unit circle                                                             
    radius = np.sqrt(1 - z * z)

    # Determine where xy fall on the sphere, given the azimuthal and polar angles                                            
    y = radius * np.sin(theta) * r
    x = radius * np.cos(theta) * r

    # Display points in a scatter plot                                                                                       
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(x, y, z)
    plt.show()


if __name__ == '__main__':
    fibonacci_sphere(int(500), 1)


#%%
N = 100

theta = np.linspace(0, 2*np.pi, 100)
r = np.linspace(0, 20**(1/6), 100) ** 6

R, THETA = np.meshgrid(r, theta)
X = R * np.cos(THETA)
Y = R * np.sin(THETA)

Z = X ** 2 + Y ** 2
plt.pcolormesh(X, Y, Z, edgecolors='k')
plt.show()

# %%
x_i1 = np.linspace(-5,5,50)
x_o1 = np.linspace(5,15,10)
# x1 = np.concatenate((x_i1, x_o1))

# x_i2 = np.linspace(-5,0,1)
x_o2 = np.linspace(-15,-5,10)
# x2 = np.concatenate((x_o2, x_i2))

Xi, Yi = np.meshgrid(x_i1, x_i1)
Ri = Xi **2 + Yi ** 2
plt.pcolormesh(Xi, Yi, Ri, edgecolors='k', cmap='plasma')
plt.show()
xo = np.concatenate((x_o2, x_o1))
Xo, Yo = np.meshgrid(xo, xo)
Ro = Xo **2 + Yo ** 2
plt.pcolormesh(Xo, Yo, Ro, edgecolors='white')
plt.show()


# %%
x_i1 = np.linspace(-5,5,20)
x_o1 = np.linspace(5,15,10)
x_o2 = np.linspace(-15,-5,10)
x = np.concatenate((x_o2, np.concatenate((x_i1,x_o1))))
X, Y = np.meshgrid(x,x)
R = X + Y
for n, row in enumerate(X):
    if Y[n,0] < -5 or Y[n,0] > 5:
        print(Y[n,0])

plt.pcolormesh(X, Y, R, edgecolors='k')
plt.show()


# %%
lim = 50
x = np.linspace(-lim, lim, 50)
X, Y, Z = np.meshgrid(x, x, x)

a = 0.5  # how far to contract
b = 0.8# how strongly to contract
c = 1 - b*np.exp(-((X/lim)**2 + (Y/lim)**2)/a**2)
x, y = c*X, c*Y
R = x **2 + y **2
norm = LogNorm(np.min(R), np.max(R))
fig = plt.figure()
ax = fig.add_subplot(projection='3d')
colors = ax.scatter(x[...,0], y[...,0], Z[...,25], edgecolors='k')
fig.colorbar(colors)
ax.set_xlabel('X')
ax.set_ylabel('y')
plt.show()



# %%
x = np.linspace(-10,10, 20)
y = np.linspace(-20,20, 25)
X, Y = np.meshgrid(x, y)
Z = np.sqrt(np.square(X) + np.square(Y))
mask = (X < 5) * (Y < 5)
mask1 = (-5 < X) * (-5 < Y)
new_Z = np.where(mask * mask1, 0, Z)

xi = np.linspace(-5,5, 30)
yi = xi
Xi, Yi = np.meshgrid(xi, yi)
Zi = np.sqrt(np.square(Xi) + (Yi))


plt.pcolormesh(X, Y, new_Z)
plt.pcolormesh(Xi, Yi, Zi)
plt.colorbar()
plt.show()

# %%
import numpy as np

wvl = np.load('G_(T,L)/geom-wvl.npy')
temps = np.load('G_(T,L)/temps.npy')
G = np.load('G_(T,L)/G-1e-05_geomwvl.npy')

norm = LogNorm(vmax=np.max(G), vmin=np.max(G)/1e6, clip=True)
L, T = np.meshgrid(wvl, temps)
plt.pcolormesh(L, T, G, norm=norm)
plt.colorbar()
# plt.yscale('log')
plt.xscale('log')
plt.xlabel('Wavelength ($\AA$)')
plt.ylabel('Temperature (K)')
plt.savefig('Figures/Hard-Xray-G(T,L).png', dpi=500)
# %%
import pandas as pd
data = pd.read_csv('/Users/tobias/Data/Extended_StellarParams.csv', dtype=None, delimiter=';', encoding='utf-8')


data_np = np.genfromtxt('/Users/tobias/Data/Extended_stellar_params.csv',delimiter=';', dtype=[('Star', 'U10'), ('Association', 'U10'), ('Age', 'f8'), ('Age_er', 'f8'),
                                                                                 ('Prot', 'f8'), ('Prot_err', 'f8'), ('Teff', 'f8'), ('Teff_err', 'f8'),
                                                                                 ('logg', 'U10'), ('vsin_i', 'U10'), ('xi', 'U10'), ('vr', 'U10')], skip_header=1, encoding='UTF-8')
names = [i for i in os.listdir(os.environ['FMPdata']) if i[0] != '.']
print(data.keys)

# %% 
a = np.zeros((3,500))
b = np.ones((500,4))
print(np.dot(a, b).shape)
print(a.shape)
a = a[:, np.newaxis]
print(a.shape)

# %%
# Calculate from L from F (ROSAT)
f1 = np.array([1.24e-11, 2.28e-11, 3.00e-11, 3.09e-12, 1.28e-11, 1.90e-11])
f2 = np.array([5.34e-12, 6.21e-12, 9.27e-12, 9.68e-13, 6.97e-13, 4.94e-12])
d = np.array([2.1e-2, 3.5e-2, 1.0e-2, 2.20e-2, 2.40e-2, 1.40e-2]) * 3.086e21
L1 = 4* np.pi * np.square(d) * f1
L2 = 4 * np.pi * np.square(d) * f2
print(L1)
print(L2)
# %%
import Grid_Operations as GO
import numpy as np
import matplotlib.pyplot as plt
grid, x= GO.create_grid(10, 90, 'linear')
segments = GO.segment_3dgrid(grid)

ax = plt.figure().add_subplot(projection='3d')

for i, (X, Y, Z) in enumerate(segments):
    print(X.shape, Y.shape, Z.shape, i)
    if i == 14:
        continue
    ax.scatter(X, Y, Z, alpha=0.1)