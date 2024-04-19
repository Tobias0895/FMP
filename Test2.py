import numpy as np
import matplotlib.pyplot as plt
from Xray_winds.src.Xray_winds.Star_class import star_model

lstar = star_model('1x-Pwand', interpolation='linear')
nstar = star_model('1x-PWand', interpolation='nearest')

npix = 91
llx = lstar.lum_x(250, npix, grid_type='segmented')
nlx = nstar.lum_x(250, 9*npix, grid_type='segmented')

print(llx/nlx)