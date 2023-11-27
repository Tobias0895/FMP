#%%
import numpy as np
from scipy import interpolate


__all__ = ['create_grid', 'rotate_grid']

def create_grid(size, resolution):
    x = np.linspace(-size/2, size/2, resolution)
    X, Y, Z = np.meshgrid(x,x,x)
    return X, Y, Z


def rotate_grid(angle:float|int, grid:tuple):
    X, Y, Z = grid
    Rotation_matrix = np.array([[np.cos(angle), np.sin(angle), 0],
                       [np.sin(angle), np.cos(angle), 0],
                       [0, 0, 1]])
    
    X_prime = X * np.cos(angle) - Y * np.sin(angle)
    Y_prime = X * np.sin(angle) + Y * np.cos(angle)
    Z_prime = Z.copy()
    return X_prime, Y_prime, Z_prime

if __name__ =='__main__':
    import matplotlib.pyplot as plt
    from scipy.interpolate import NearestNDInterpolator
    import os 
    from starwinds_readplt.dataset import Dataset

    data_loc = os.environ['FMPdata']
    ds = Dataset.from_file(data_loc + '/3d__var_3_n00060000.plt')
    ds_data = np.stack([ds(name) for name in ds.variables], axis=-1)
    ds_points = ds_data[...,:3]

    nearest = NearestNDInterpolator(ds_points, ds_data)
    angels = np.linspace(0, 350, 50)
    
    X, Y, Z = create_grid(20, 100)
    data = nearest(X, Y, Z)
    for angle in angels:
        X_rot, Y_rot, Z_rot = rotate_grid(np.deg2rad(angle), (X, Y, Z))
        data_prime = nearest(X_rot, Y_rot, Z_rot)
        fig = plt.figure()
        ax = fig.add_subplot(121, projection='3d')
        surf = ax.plot_surface(X_rot[0,...], Z_rot[0,...], data_prime[:, 10, :,ds.variables.index('te [K]')])
        # plt.colorbar()
        plt.show()

