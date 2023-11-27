from moviepy.editor import VideoClip
from moviepy.video.io.bindings import mplfig_to_npimage
import numpy as np
import matplotlib.pyplot as plt
import Test_site as ts

__all__ = ['make_frame']

def make_frame(t, **kwargs):
    from Test_site import nearest_interpolater
    fig = plt.figure()
    ax = fig.add_subplot()
    ax.clear()

    img, mesh = ts.integrate_along_line_of_sight((20,50), direction='+x', stellar_radius=1, interpolater=nearest_interpolater,
                                                    angle=np.deg2rad(t), image_radius=20, pixel_count=100)
    ax.pcolormesh(mesh[0], mesh[1], img, log='norm')

    return mplfig_to_npimage(fig)


if __name__ == '__main__':
    movie = VideoClip(make_frame, duration=300)
    movie.ipython_display(fps = 20, loop=True, auto_play=True)