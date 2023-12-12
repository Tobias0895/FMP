# %%
from matplotlib.gridspec import GridSpec
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from mpl_toolkits.axes_grid1 import make_axes_locatable

# %% 
# Make a plot of G(T,L)
G = np.load('G_(T,L)/G-1e-05.npy')
wvl = np.load('G_(T,L)/L.npy')
temps = np.load('G_(T,L)/T.npy')

vmax = np.max(G)
norm = LogNorm(vmin=vmax/1e6, vmax=vmax, clip=True)

fig = plt.figure(figsize=(10,10), layout='constrained')
gs = GridSpec(3,3, hspace=0, wspace=0, width_ratios=[1,1,3], height_ratios=[3,1,1], )



ax1 = fig.add_subplot(gs[2,1:])
ax1.plot(wvl[100, 100:], G[100, 100:]/np.max(G[100, 100:]), label=f'{round(temps[100,0])} K')
ax1.plot(wvl[150, 100:], G[150, 100:]/np.max(G[100, 100:]), label=f'{round(temps[150,0])} K')
ax1.set_yscale('log')
ax1.set_xlabel('Wavelength [$\AA$]')
ax1.legend()


ax2 = fig.add_subplot(gs[0:2,0] )
ax2.plot(G[:,1500]/np.max(G[:,1500]), temps[:,1500], label=f'$\lambda$ = {wvl[0,1500]}')
ax2.plot(G[:,500]/np.max(G[:,1500]), temps[:,500], label=f'$\lambda$ = {wvl[0,500]}')
ax2.legend()
ax2.set_yscale('log')
ax2.set_ylabel('T [K]')

ax0 = fig.add_subplot(gs[0:2,1:], sharex=ax1, sharey=ax2)
colors = ax0.pcolormesh(wvl, temps, G, norm=norm)
ax0.set_xticklabels
cax = inset_axes(ax0, width='6%', height='100%', borderpad=-4, loc='right')
fig.colorbar(colors, extend='min', cax=cax, label='erg cm$^{3}$ s$^{-1}$ sr$^{-1}$ $\AA^{-1}$')
plt.setp(ax0.get_xticklabels(), visible=False)
plt.setp(ax0.get_yticklabels(), visible=False)

plt.show()
