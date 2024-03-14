from Calculate_flux import projection_2d
import Star_class
import numpy as np
import matplotlib.pyplot as plt
import load_data
from scipy.optimize import curve_fit

name = '1x-HHLeo'
star = Star_class.star_model(name)
interpolator, var_list, star_params, model_params, ds_points, ds_data = load_data.import_data(name, full_output=True)
# Create a grid
size = 20 # stellar radii
x = np.linspace(-size, size, 100) * star_params["RadiusStar"]
X, Y, Z = np.meshgrid(x,x,x)
XX, YY = np.meshgrid(x,x)
Z = np.full_like(XX, 1) * star_params['RadiusStar']

R = np.array([XX, YY, Z]).T
interpolated_data = interpolator(XX, YY, Z)
B_vector = np.asarray([interpolated_data[...,var_list.index('B_x [Gauss]')],
            interpolated_data[...,var_list.index('B_y [Gauss]')]]).T


def field_function2(r_vector, m_vector):
    mu_0 =  1#1.25663706212e-6
    r_size = np.linalg.norm(r_vector, axis=-1)
    r_size = r_size[...,np.newaxis]
    r_hat = r_vector / r_size
    a = 3 * (np.dot(r_hat, m_vector))
    b = (r_hat - m_vector.T)
    c = (a * b) / (np.power(r_size, 3))
    return mu_0 / (4 * np.pi) * c


m_vector = np.array([.1,-0.2,0])
m_vector = m_vector[...,np.newaxis]
field2 = field_function2(R, m_vector)

fig, ax = plt.subplots(1,2, figsize=(16,9))
ax[0].streamplot(XX/star_params['RadiusStar'], YY/star_params['RadiusStar'], field2[...,0], field2[...,1])
ax[1].streamplot(XX/star_params['RadiusStar'], YY/star_params['RadiusStar'], B_vector[...,0], B_vector[...,1])
plt.show()


