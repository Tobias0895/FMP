from matplotlib.figure import figaspect
from matplotlib.transforms import Bbox
from moviepy.editor import VideoClip
from moviepy.video.io.bindings import mplfig_to_npimage
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
import Calculate_flux as cf

__all__ = ['make_lightcurve_frame', 'make_rotation_frame']


def make_lightcurve_frame(t):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.clear()
    current_angle = t * 36
    wvl, spectrum = cf.create_spectra((1,180), 15, disable_tqdm=True, stellar_radius=1, direction='+x',
                                        interpolator=interpolator, var_list=v, angle=np.deg2rad(current_angle))
    ax.plot(wvl, spectrum)
    ax.set_xlabel('Angle (deg)')
    ax.set_ylabel('Flux')

    return mplfig_to_npimage(fig)
    


def make_rotation_lightcurve_frame(t, **kwargs):
    from matplotlib.colors import LogNorm
    fig = plt.figure()
    angle = t * 36
    fig, ax = plt.subplots(2, 1, figsize=(5,7), gridspec_kw={'height_ratios':[3,1]})
    ax[0].clear()
    ax[1].clear()
    img, mesh = cf.projection_2d(band, direction='+x', stellar_radius=params['RadiusStar'], interpolator=interpolator,
                                                    angle=(np.deg2rad(angle),0), image_radius=20, pixel_count=100, var_list=v, use_simple_g=False)
    vmax = np.nanmax(img)
    norm = LogNorm(vmax=vmax, vmin=vmax/1e6 , clip=True)
    color = ax[0].pcolormesh(mesh[0], mesh[1], img, norm=norm)
    ax[0].set_title('Angle = {:.1f}'.format(angle))
    ax[0].axis('equal')
    fig.colorbar(color)
    ax[0].set_ylabel('Z [R]')
    ax[0].set_xlabel('d [R]')

    angle_rad = np.deg2rad(angle)
    light_curve =np.load(f'light_curves/Light_curve_{name}.npy')
    angles, lum = (light_curve[0], light_curve[1])
    lum_at_angle = np.interp(angle_rad, angles, lum)
    ax[1].plot(angles, lum)
    ax[1].scatter(angle_rad, lum_at_angle)
    ax[1].set_xlabel('Angle (rad)')
    ax[1].set_ylabel("L_x erg s^-1")

    return mplfig_to_npimage(fig)


def make_rotation_frame(t):
    from matplotlib.colors import LogNorm
    fig = plt.figure()
    angle = t * 36
    fig, ax = plt.subplots(1, 1)
    ax.clear()
    img, mesh = cf.projection_2d(band, direction='+x', stellar_radius=params['RadiusStar'], interpolator=interpolator,
                                                    angle=(np.deg2rad(angle), 0), image_radius=20, pixel_count=100, var_list=v, use_simple_g=False)
    vmax = np.nanmax(img)
    norm = LogNorm(vmax=vmax, vmin=vmax/1e6 , clip=True)
    color = ax.pcolormesh(mesh[0], mesh[1], img, norm=norm, shading='gouraud')
    ax.set_title('Angle = {:.1f}'.format(angle))
    ax.axis('equal')
    fig.colorbar(color)
    ax.set_ylabel('Z [R]')
    ax.set_xlabel('d [R]')

    return mplfig_to_npimage(fig)

if __name__ == '__main__':
    import load_data
    name = '1x-Mel25-005'
    interpolator, v, params = load_data.import_data(name) # type: ignore
    band = (10, 30)
    # spectrum_movie = VideoClip(make_lightcurve_frame, duration=10)
    # spectrum_movie.write_gif('1x-Mel25-005_spectrum.gif', fps=24)

    rotation_movie = VideoClip(make_rotation_lightcurve_frame, duration=10)
    rotation_movie.write_gif(f'Movies/rotation_curve_{name}.gif', fps=24)