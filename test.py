import Grid_Operations as GO
import numpy as np
import matplotlib.pyplot as plt
import Star_class
import Calculate_flux as cf

star = Star_class.star_model('1x-PWAnd')
print(star.lum_x(50, 100))

npix = 100
grid, x= GO.create_grid(50 * star.params['RadiusStar'], npix, 'linear')

r_star = 5
segments = GO.up_center_res(grid, 3)
ax = plt.figure().add_subplot()

total_X = []
for i, segment in enumerate(segments):
    X, Y, Z = segment
    X_rot, Y_rot, Z_rot = GO.rotate_grid(0, 0, segment)

    star_mask = X ** 2 + Y ** 2 + Z ** 2 <= r_star ** 2
    interpolated_data = star.interpolator(X_rot, Y_rot, Z_rot)
    integrand = np.square(interpolated_data[...,star.var_list.index('Rho [g/cm^3]')] / 1.67e-24) * cf.G(interpolated_data[...,star.var_list.index('te [K]')], (0.1, 180))
    masked_integrand = np.where(star_mask==False, integrand, 0)

    area_of_star =  np.square(r_star) * np.pi
    solid_angle_array = area_of_star / (X_rot ** 2 + Y_rot ** 2 + Z_rot ** 2)
    fraction_light = 1 - (solid_angle_array / (4*np.pi))
    masked_integrand *= fraction_light

    projection = np.trapz(masked_integrand, X, axis=1)

    one_d = np.trapz(projection, Z[:, 0 , :], axis=-1)
    tot = np.trapz(one_d, Z[0, 0, :])
    total_X.append(tot)

print(np.sum(total_X))
plt.show()