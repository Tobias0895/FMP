import numpy as np
import matplotlib.pyplot as plt

def read_data(file, path='Data/', print_header=False):

    from starwinds_readplt.dataset import Dataset
    ds = Dataset.from_file(path + file)
    if print_header:
        print(ds)
    return ds

def generate_spectra(wvl, temperatures, densities, Ion='h_1s', plot=False):
    import ChiantiPy.core as ch
    print(densities[25,25], temperatures[25,25])
    densities = np.ravel(densities)
    temperatures = np.ravel(temperatures)
    number_densities = densities * 6.02e23
    c = ch.continuum(Ion, temperature=temperatures, em=1e40)
    c.freeFree(wvl)

    h_1 = ch.ion(Ion, temperature=temperatures, eDensity=number_densities, em=1e+27)
    h_1.spectrum(wvl)
    Total_spectrum = c.FreeFree['intensity'] #+ c.FreeBound['intensity']
    return Total_spectrum



if __name__ == '__main__':
    from scipy.interpolate import NearestNDInterpolator
    ds = read_data('3d__var_3_n00060000.plt', print_header=False)
    ds_points = np.stack([ds(name) for name in ds.variables[:3]], axis=-1)
    ds_data = np.stack([ds(name) for name in ds.variables], axis=-1)
    interpolator = NearestNDInterpolator(ds_points, ds_data)

    resolution = 60
    x = np.linspace(-30, 30, resolution)  # Units are solar/stellar radii
    y = np.linspace(-30, 30, resolution +  1)  # Use different number of elements in x and y direction to avoid confusing them later on.
    X, Y = np.meshgrid(x, y)
    Z = np.full_like(X, 0)

    interpolated_data = interpolator(X, Y, Z)
    wavelength = np.linspace(0.1, 100, 1000) # 0.1 - 100 Xray range in angstrom
    spectra = generate_spectra(wavelength, interpolated_data[...,ds.variables.index('ti [K]')],
                                interpolated_data[...,ds.variables.index('Rho [g/cm^3]')])
    spectra = spectra.reshape(interpolated_data.shape[0], interpolated_data.shape[1], len(wavelength))

    plt.plot(wavelength, spectra[25,25,:], label='middle pixel')

    plt.legend()
    plt.xlabel("Wavelength ($\AA$)")
    plt.ylabel("erg cm$^{-2}$ s$^{-1}$ $\AA^{-1}$ Sr$^{-1}$")
    # plt.ylim(0, 1.4e-30)
    plt.show()