from Xray_winds.src.Xray_winds import Grid_Operations as GO
from Xray_winds.src.Xray_winds import Star_class
import matplotlib.pyplot as plt 
import numpy as np

star = Star_class.star_model('1x-PWAnd')
print(star.lum_x(250, 121, grid_type='segmented', nseg=1))


npix = 121
print((npix - 1) % 3)
x = np.linspace(-10,10,npix)
grid = np.meshgrid(x,x,x, indexing='ij')

segments = GO.up_center_res(grid)
# middle = segments.pop(-1)
# segments_inner = GO.up_center_res(middle)
# segments += segments_inner
seg_x_tot = []
segment_sum = []
ax = plt.figure().add_subplot(projection='3d')
for i, segment in enumerate(segments):
    X, Y, Z = segment
    masked_integrand = np.full_like(X, 1)
    projection = np.trapz(masked_integrand, X, axis=0)
    # print(projection.max())
    one_d = np.trapz(projection, Y[0,:,:], axis=0)
    seg_tot = np.trapz(one_d, Z[0, 0, :])
    segment_sum.append(np.sum(masked_integrand))
    seg_x_tot.append(seg_tot)

print(f"intgral: {np.sum(seg_x_tot)}")
lin_grid, x= GO.create_grid(10, 100, 'linear')
Xl, Yl, Zl = lin_grid
full_cube = np.full_like(Xl, 1)
two_d = np.trapz(full_cube, Xl, axis=0)
one_d = np.trapz(two_d, Yl[0,:,:], axis=0)
tot = np.trapz(one_d, Zl[0, 0, :])
print(f"integral: {tot}")
print(f"segmented gives {100 * np.sum(seg_x_tot) / tot:.2f}% from linear")
