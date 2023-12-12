#%%
import numpy as np
import ChiantiPy.core as ch
import matplotlib.pyplot as plt
from ChiantiPy.tools import filters
from matplotlib.colors import LogNorm
import random

time = np.linspace(0, 10, 40)
for i in time:
    max_angle = 360
    current_angle = i * 36
    print(int(current_angle/max_angle * 40))
    angles = np.linspace(0, current_angle, )


# %%
import os
data_loc = os.environ['FMPdata']
model_name = '1x-PWAnd'
star_file = open(os.path.join(data_loc, model_name + '/STAR.in'), 'r')
param_lines = star_file.readlines()
star_params = {}
for n, line in enumerate(param_lines[-3:]):
    split = line.split('  ')
    value = split[0]
    name = split[-1].split(' ')
    if name[0].isalpha() == False:
        name = split[-1].split(' ')[1]
    else:
        name = name[0]
    star_params[str(name)] = value
print(star_params['RotationPeriodStar'])


# %%
def count(*args):
    for i in args:
        print(i)


count(*[1,3,1,4,1,1,1,1])