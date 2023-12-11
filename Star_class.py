import matplotlib
import load_data
import matplotlib.pyplot as plt
import numpy as np
import Calculate_flux
import matplotlib as mpl

class star_model():

    def __init__(self, name, interpolation='nearest'):
        self.name = name
        data = load_data.import_data(self.name, interpolate=interpolation)
        self.interpolator = data[0] 
        self.var_list = data[1]
        self.params = data[2]

    def raw_data(self):
        return load_data.read_model(self.name)
    
    def pop_plot_star(self, theta:float|int, phi, wavelength_range:tuple, save='', **grid_kw):
        from matplotlib.colors import LogNorm
        lum, mesh = Calculate_flux.projection_2d(wavelength_range, self.params['RadiusStar'], self.interpolator, self.var_list, angle=(theta, phi), **grid_kw)
        vmax = np.nanmax(lum)
        norm = LogNorm(vmax=vmax, vmin=vmax/1e6)
        fig = plt.figure()
        ax = fig.add_subplot(111)

        quadmesh = ax.pcolormesh(mesh[0], mesh[1], lum, norm=norm, shading='gouraud')
        
        fig.colorbar(quadmesh)
        ax.set_xlabel('r (cm)')
        ax.set_ylabel('r (cm)')
        
        if save:
            plt.savefig(save, dpi=500)
        plt.show()
        return fig

    def spectrum(self, wvl_range:tuple, spec_res:int|float, **spec_kwargs) -> tuple:
        wvl, spectrum = Calculate_flux.create_spectra(wvl_range, spec_res, stellar_radius=self.params['RadiusStar'], interpolator=self.interpolator,
                                                    var_list=self.var_list, **spec_kwargs)
        return wvl, spectrum
    
    def spectrum_pop_plot(self, wvl_range:tuple, spec_res:int|float, **spec_kwargs):
        wvl, spec = self.spectrum(wvl_range, spec_res, **spec_kwargs)
        spec_fig = plt.figure()
        ax = spec_fig.add_subplot(111)
        ax.plot(wvl, spec)
        ax.set_yscale('log')
        ax.set_xlabel('Wavelength ($\AA$)')
        ax.set_ylabel('erg s$^{-1}$ $\AA^{-1}$')
        plt.show()



    def light_curve(self, rotation_direction:str, wvl_range=(0.1, 180)):
        from tqdm import tqdm
        angles = np.linspace(0, 2*np.pi, 50) # in radians
        fluxes = np.zeros(len(angles))
        if rotation_direction.lower() == 'equator':
            for i, ang in enumerate(tqdm(angles)):
                fluxes[i] = Calculate_flux.total_lum_wvl_bin(wvl_range, stellar_radius=self.params['RadiusStar'], interpolator=self.interpolator,
                                                              var_list=self.var_list, angle=(ang, 0.))
        elif rotation_direction.lower() == 'poles':
            for i, ang in enumerate(tqdm(angles)):
                fluxes[i] = Calculate_flux.total_lum_wvl_bin(wvl_range, stellar_radius=self.params['RadiusStar'], interpolator=self.interpolator, 
                                                             var_list=self.var_list, angle=(0., ang))
        return angles, fluxes
    
    def lum_x(self, wvl_bin=(0.1,180)):
        return Calculate_flux.total_lum_wvl_bin(wvl_bin, self.params['RadiusStar'], self.interpolator, var_list=self.var_list)

if __name__ == "__main__":
    sun = star_model('sun')
    print(np.log10(sun.lum_x()/(3.828e33)))
    