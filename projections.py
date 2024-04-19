from matplotlib.axis import YTick
from Xray_winds.src.Xray_winds import Star_class
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import numpy as np
import os, sys
from tqdm import tqdm

names = [i for i in os.listdir(os.environ['FMPdata']) if i[0] != '.']
fig = plt.figure(figsize=(15,6), layout='constrained')
gs = GridSpec(2, 6, width_ratios=[1,1,1,1,1,0.1], wspace=0.1, hspace=0.2)
for i in tqdm(range(5)):
    ax = fig.add_subplot(gs[0, i])
    ax2 = fig.add_subplot(gs[1, i])
    star = Star_class.star_model(names[i+1])
    star2 = Star_class.star_model(names[-1-i])
    outermesh = star.projection_figure(0,0,[5,120], ax=ax, grid_type='segmented')
    outermesh2 = star2.projection_figure(0,0,[5,120], ax=ax2, grid_type='segmented')
    ax.set_yticklabels([]); ax.set_xticklabels([])
    ax.set_title(names[i+1])
    ax2.set_title(names[-1-i])
    if i > 0:
        plt.setp([ax, ax2], yticklabels=[])
        plt.setp([ax, ax2], xticklabels=[])

cax1 = fig.add_subplot(gs[:, -1])
fig.colorbar(outermesh, cax=cax1)
fig.savefig('Figures/Paper/Projection_comp.pdf', dpi=500, bbox_inches='tight')
plt.show()

