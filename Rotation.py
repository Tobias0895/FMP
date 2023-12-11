#%%
import numpy as np
from scipy import interpolate


__all__ = ['create_grid', 'rotate_grid']

def create_grid(size, resolution:int, type:str):
    if type.lower() == 'linear':
        x = np.linspace(-size/2, size/2, resolution)
        X, Y, Z = np.meshgrid(x,x,x)
        return (X, Y, Z)
    
    elif type.lower() == 'log':
        # create the positive half
        x_p = np.logspace(0, np.log10(size/2), resolution)
        x_n = -1 * x_p
        x = np.append(np.flip(x_n), x_p)
        X, Y, Z = np.meshgrid(x,x,x)
        return (X, Y, Z)
    

def rotate_grid(theta:float|int, phi:float|int, grid:tuple):
    X, Y, Z = grid

    # First rotate around Z-axis (Theta)
    X_prime = X * np.cos(theta) - Y * np.sin(theta)
    Y_prime = X * np.sin(theta) + Y * np.cos(theta)
    Z_prime = Z.copy()

    # Then rotate around Y (phi)
    X_prime_prime = X_prime * np.cos(phi) + Z_prime * np.sin(phi)
    Y_prime_prime = Y_prime.copy()
    Z_prime_prime = X_prime * - np.sin(phi) + Z_prime * np.cos(phi)

    return X_prime_prime, Y_prime_prime, Z_prime_prime

if __name__ == '__main__':
    import matplotlib.pyplot as plt
    from scipy.interpolate import NearestNDInterpolator
    import os 
    from starwinds_readplt.dataset import Dataset
    import Calculate_flux as cf


    nearest, var_list = cf.import_data('sun', 1)
    
    X, Y, Z = create_grid(20, 60, 'Log')

    X_rot, Y_rot, Z_rot = rotate_grid(np.deg2rad(90), (X, Y, Z))

    data_prime = nearest(X_rot, Y_rot, Z_rot)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    surf = ax.pcolormesh(X_rot[0,...], Z_rot[0,...], data_prime[10, :, :,var_list.index('te [K]')])
    # plt.colorbar()
    plt.show()

