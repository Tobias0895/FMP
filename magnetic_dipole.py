import numpy as np
import load_data
import Star_class
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit, minimize, basinhopping

name = '1x-PWAnd'

interpolator, var_list, star_params, model_params, ds_points, ds_data = load_data.import_data(name, full_output=True)
# Create a grid
mask = np.sum(np.square(ds_points/star_params['RadiusStar']), axis=-1) <= (3)

R_points = np.array([ds_points[...,0][mask], ds_points[...,1][mask], ds_points[...,2][mask]]).T
B_points = np.array([ds_data[...,var_list.index('B_x [Gauss]')][mask], 
            ds_data[...,var_list.index('B_y [Gauss]')][mask],
            ds_data[...,var_list.index('B_z [Gauss]')][mask]]).T

print(R_points.shape)
x = np.linspace(4, 4, 10) #* star_params['RadiusStar'] #* 1e-2 # in meters
X, Y, Z = np.meshgrid(x,x,x)
interpolated_data = interpolator(X, Y, Z)
pos = np.array([X, Y, Z], dtype=np.float64).T
true_field = np.asarray([interpolated_data[...,var_list.index('B_x [Gauss]')], 
            interpolated_data[...,var_list.index('B_y [Gauss]')],
            interpolated_data[...,var_list.index('B_z [Gauss]')]], dtype=np.float64).T
true_field /= np.max(true_field) # normalize the field 

def dipfield(r_vector, m_vector):
    mu_0 = 1 #1.25663706212e-6
    r_size = np.linalg.norm(r_vector, axis=-1)
    r_size = r_size[...,np.newaxis]
    r_hat = r_vector / r_size
    a = 3 * (np.dot(r_hat, m_vector))
    b = (a * r_hat - m_vector.T)
    c = b / (np.power(r_size, 3))
    return mu_0 / (4 * np.pi) * c

def error(m_vector):
    m_vector = m_vector[...,np.newaxis]
    field = dipfield(R_points, m_vector) * 1e5
    field /= np.max(field) # normalize the field
    return np.sum(np.square(field - B_points))

def two_component_error(m_vector):
    mx, my = m_vector
    mz = np.sqrt(1 - mx**2 - my**2)
    m_vector = np.asarray([mx, my, mz])[...,np.newaxis]
    field = dipfield(pos, m_vector) * 1e5 # in Teslas
    field /= np.max(field) # normalize the field
    return np.sum(np.square(field - true_field))

def magnetic_inclination(m_vector):
    """Calculates the magnetic inclination. In all the models the rotation axis is alligned with the Z-axis
    so the inner product of M * Z = |M| * |Z| cos(theta) -> theta = arccos((M * Z) / (|M|*|Z|)).
    since |Z| = 1, theta = arccos((m_vector * Z_normal) / |M|)

    Args:
        m_vector (list:ndarray): The dipole moment that discribes the first order approximation of the field
    """
    z_normal = np.array([0,0,1], dtype=np.float64)
    return np.rad2deg(np.arccos(np.dot(m_vector, z_normal) / np.linalg.norm(m_vector)))
    


m_vector = np.array([-2e3, 1e2,-1e4], dtype=np.float64)
res1 = minimize(error, m_vector, tol=1e-8) #, method='CG', options={'disp':True, 'gtol':1e-10})
res2 = minimize(error, m_vector, method='CG', options={'disp':True, 'gtol':1e-10})
res3 = minimize(two_component_error, [-2e3, 1e4], tol=1e-8)
x3 = np.array([res3['x'][0], res3['x'][1], np.sqrt(1 - res3['x'][0]**2 - res3['x'][1]**2)])
print(res3)
print(res1)
print(magnetic_inclination(res2['x']), 180 - magnetic_inclination(res2['x']))
print(magnetic_inclination(x3), 180 - magnetic_inclination(x3))   