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
        X, Y = np.meshgrid(x,x)
        Z = X ** 2 + Y ** 2
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
    import load_data


    nearest, var_list, params = load_data.import_data('sun')
    
    X, Y, Z = create_grid(20 * params['RadiusStar'], 60, 'Log')

    print('grid made')
    fig = plt.figure()
    ax = fig.add_subplot(111)
    surf = ax.pcolormesh(X, Y, Z)
    # plt.colorbar()p[]
    plt.show()

