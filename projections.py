# %%
from matplotlib.axis import YTick
import matplotlib.style
from Xray_winds.src.Xray_winds import Star_class
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import numpy as np
import os, sys
from tqdm import tqdm
import matplotlib 


matplotlib.style.use('mnras.mplstyle')
# %%
mnras_single_width_in = 3.3209
mnras_double_width_in = 6.9738
names = [i for i in os.listdir(os.environ['FMPdata']) if i[0] != '.']
print(names)
names = ['Sun']
fig = plt.figure(figsize=(mnras_double_width_in, 3/4 * mnras_single_width_in), layout='constrained')
gs = GridSpec(2, 6, width_ratios=[1,1,1,1,1,0.1], wspace=0.1, hspace=0.2)
for i in tqdm(range(5)):
    ax = fig.add_subplot(gs[0, i])
    ax2 = fig.add_subplot(gs[1, i])
    star = Star_class.star_model(names[i+1])
    # star2 = Star_class.star_model(names[-1-i])
    outermesh = star.projection_figure(0,0,[0.,180], ax=ax, grid_type='segmented', vmax=1e5, image_radius=5)
    # outermesh2 = star2.projection_figure(0,0,[0,180], ax=ax2, grid_type='segmented', vmax=1e5, image_radius=5)
    ax.set_yticklabels([]); ax.set_xticklabels([])
    # ax.set_title(star.name)
    # ax2.set_title(star2.name)

    ax.text(-4.5, 3.5, f'{star.name} \n P$_{{rot}}$ = {star.params["RotationPeriodStar"]}', color='white')
    # ax2.text(-4.5,3.5, f'{star2.name} \n P$_{{rot}}$ = {star2.params["RotationPeriodStar"]}', color='white')
    if i > 0:
        plt.setp([ax, ax2], yticklabels=[])
        plt.setp([ax, ax2], xticklabels=[])
    break
cax1 = fig.add_subplot(gs[:, -1])
ax = fig.axes
ax[1].set_xlabel('Y [R$_{*}$]')
ax[1].set_ylabel('Z [R$_{*}$]')

fig.colorbar(outermesh, cax=cax1)
fig.savefig('Figures/Size_test_fig.pdf', dpi=500, bbox_inches='tight')
plt.show()


# %%

# %%
4/5
# %%
