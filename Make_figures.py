# %%
from matplotlib.gridspec import GridSpec
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import matplotlib.style
import numpy as np
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib
matplotlib.rcParams.update(matplotlib.rcParamsDefault)
matplotlib.style.use('mnras.mplstyle')

# %% 
# ------------Make a plot of G(T,L)-----------------

mnras_single_width_in = 3.3209
G = np.load('G_(T,L)/G-1e-05.npy')
wvl = np.load('G_(T,L)/L.npy')
temps = np.load('G_(T,L)/T.npy')

vmax = np.max(G)
norm = LogNorm(vmin=vmax/1e6, vmax=vmax, clip=True)

fig = plt.figure(figsize=(mnras_single_width_in, mnras_single_width_in), layout='constrained')
gs = GridSpec(3,3, hspace=0, wspace=0, width_ratios=[1,1,3], height_ratios=[3,1,1])



ax1 = fig.add_subplot(gs[2,1:])
ax1.plot(wvl[100, 100:], G[100, 100:]/np.max(G[100, 100:]), label=f'{round(temps[100,0]):.1e} K')
ax1.plot(wvl[150, 100:], G[150, 100:]/np.max(G[100, 100:]), label=f'{round(temps[150,0]):.1e} K')
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
colors = ax0.pcolormesh(wvl, temps, G, norm=norm, cmap='plasma', rasterized=True)
# tax0 = ax2.twinx(ax2.xaxis.convert_units
ax0.set_xticklabels
cax = inset_axes(ax0, width='6%', height='100%', borderpad=-3.5, loc='right')
fig.colorbar(colors, extend='min', cax=cax, label='erg cm$^{3}$ s$^{-1}$ sr$^{-1}$ $\AA^{-1}$')
plt.setp(ax0.get_xticklabels(), visible=False)
plt.setp(ax0.get_yticklabels(), visible=False)
# plt.savefig('Figures/Paper/Contribution_function_rasterized_B.pdf', dpi=500)
plt.show()

# %%

x = np.linspace(-20, 20, 100)
X, Y, Z = np.meshgrid(x, x, x)

star_mask = X ** 2 + Y ** 2 + Z ** 2 <= 10

fig = plt.figure()
ax = fig.add_subplot(projection='3d')
ax.voxels(star_mask)

fig.savefig('Figures/Paper/Star_mask.png', dpi=100)
plt.show()

# %%
# ------------ Make L_x vs rotation parameter plot -------
import pickle
with open('stars_data.pkl', 'rb') as f:
    stars = pickle.load(f)

# %% 
# -- Plotting --

from pyvo import registry
import matplotlib.colors as mcolors
CATALOGUE = "J/ApJ/743/48"
catalogue_ivoid = f"ivo://CDS.VizieR/{CATALOGUE}"
voresource = registry.search(ivoid=catalogue_ivoid)[0]
tables = voresource.get_tables()
table_name = "J/ApJ/743/48/catalog"

tap_service = voresource.get_service("tap")
tap_records = voresource.get_service("tap").run_sync(
    f'select  * from "{table_name}"',
)
 
dfw = tap_records.table.to_pandas()
print("keys: ", ", ".join(dfw.keys()))
sorted_R = np.sort(dfw['R*'])
low_R = sorted_R[:len(sorted_R)//3]
high_R = sorted_R[2*len(sorted_R)//3:]
mid_R = sorted_R[len(sorted_R)//3:2*len(sorted_R)// 3]
associations = np.unique(stars['association'])
fig_pl, ax = plt.subplots(1,2,figsize=(16,8), sharey=True)
fig_kpr, ax2 = plt.subplots(3,1, figsize=(4,12), sharex=True)
fig_fit, ax3 = plt.subplots(1,1, figsize=(8,8))

ax[0].scatter(np.log10(dfw['Prot']), dfw['Lx/bol'], color='grey', s=2)
ax[1].scatter(np.log10(1.86e-3 * dfw['R*']**-4 *dfw['Prot'] **-2), dfw['Lx/bol'], c='gray', s=2, label='Wright+2011')
ax3.scatter(np.log10(1.86e-3 * dfw['R*']**-4 *dfw['Prot'] **-2), dfw['Lx/bol'], c='gray', s=2, label='Wright+2011')

for Age in np.unique(np.sort(stars['Age'])):
    star_current_age= np.where(stars['Age'] == Age)
    for star in star_current_age:
        idx = np.where(stars['name'] == star)
        ax[0].scatter(np.log10(stars['Prot'][idx]), np.log10((10 ** 4) * (stars['ROSATx'][idx])/stars['Lbol'][idx]), marker='*',s=100, label=stars['association'][idx],
                    c=list(mcolors.TABLEAU_COLORS)[np.where(associations == stars['association'][idx])[0][0]])

        ax2[2].scatter(stars['kpr'][idx], np.log10(stars['EUV'][idx]/stars['Lbol'][idx]), marker='*', label=stars['association'][idx],
                    c=list(mcolors.TABLEAU_COLORS)[np.where(associations == stars['association'][idx])[0][0]])

        ax2[1].scatter(stars['kpr'][idx], np.log10(stars['ROSATx'][idx]/stars['Lbol'][idx]), marker='*', label=stars['association'][idx],
                    c=list(mcolors.TABLEAU_COLORS)[np.where(associations == stars['association'][idx])[0][0]])
        
        ax2[0].scatter(stars['kpr'][idx], np.log10(stars['HARDx'][idx]/stars['Lbol'][idx]), marker='*', label=stars['association'][idx],
                    c=list(mcolors.TABLEAU_COLORS)[np.where(associations == stars['association'][idx])[0][0]])
        
        ax[1].scatter(stars['kpr'][idx], np.log10((10 ** 4) * stars['ROSATx'][idx]/stars['Lbol'][idx]), marker='*', s=100, label=stars['association'][idx][0],
                    c=list(mcolors.TABLEAU_COLORS)[np.where(associations == stars['association'][idx])[0][0]])
        ax3.scatter(stars['kpr'][idx], np.log10((10 ** 4) * stars['ROSATx'][idx]/stars['Lbol'][idx]), marker='*', s=100, label=stars['association'][idx][0],
                    c=list(mcolors.TABLEAU_COLORS)[np.where(associations == stars['association'][idx])[0][0]])

ax[0].set_ylabel('Log L$_x$/L$_{bol}$ (erg s$^{-1}$)')
ax[0].set_xlabel('Log P$_{rot}$ (Days)')
fig_pl.subplots_adjust(wspace=0.03)

ax2[0].invert_xaxis() ; ax2[1].invert_xaxis(); ax2[2].invert_xaxis()
# ax2[0].set_ylabel('Log L$_x$ / L$_{bol}$ (erg s$^{-1}$)')
# ax2[0].axvline(-3.14, ls='--', color='grey') ; ax2[1].axvline(-3.14, ls='--', color='grey') ; ax2[2].axvline(-3.14, ls='--', color='grey')
ax2[2].set_xlabel('Log k P$^{-2}$ R$^{-4}$') # ; ax2[1].set_xlabel('k P$^{-2}$ R$^{-4}$') ; ax2[2].set_xlabel('k P$^{-2}$ R$^{-4}$')
ax2[0].text(-2,-20.5,'Hard Xrays [0.1-5] $\AA$') ; ax2[1].text(-2,-9.7, 'ROSAT [5-120] $\AA$') ; ax2[2].text(-2, -9.4, 'EUV [120-180] $\AA$')
fig_kpr.text(-0.08, 0.5, 'Log L$_x$ / L$_{bol}$ (erg s$^{-1}$)', va='center', rotation='vertical')
fig_kpr.subplots_adjust(hspace=0.05)

ax[1].invert_xaxis()
# ax[1].set_ylabel('Log L$_x$ / L$_{bol}$ (erg s$^{-1}$)')
ax[1].set_xlabel('Log k P$^{-2}$ R$^{-4}$')
# ax[1].set_title('ROSAT [5-120] $\AA$')
ax3.invert_xaxis()
kpr_axs = fig_pl.gca()
handles, labels = kpr_axs.get_legend_handles_labels()
by_label = dict(zip(labels, handles))
ax[1].legend(by_label.values(), by_label.keys())
ax2[1].legend(by_label.values(), by_label.keys(), prop={'size': 9, }, loc=6)
# fig_lit.savefig('Figures/Paper/Segmented_KPR_fig.pdf', dpi=100)
fig_kpr.savefig('Figures/Paper/KPR_in_bands.pdf', dpi=100, bbox_inches='tight')
fig_pl.savefig('Figures/Paper/Pl_Kpr_together.pdf', dpi=100, bbox_inches='tight')


# %% 
# ------------- color-color plot ------------

fig, ax = plt.subplots(1,1, figsize=(8,8))
ax.scatter(np.log10(stars['HARDx']) - np.log10(stars['ROSATx']), np.log10(stars['ROSATx']) - np.log10(stars['EUV']),
           c='k')
ax.set_xlabel("Hard - Rosat")
ax.set_ylabel("Rosat - EUV")
fig.savefig('Figures/Color_plot.png', bbox_inches='tight')
plt.show()
# %% 
# ------------- Plot stars ---------------------
import os
import Star_class
names = [i for i in os.listdir(os.environ['FMPdata']) if i[0] != '.']

for name in names:
    ang, lum = np.load(f'light_curves/Light_curve_20R_{name}.npy')
    min_lum_angle = ang[np.argmin(lum)]
    star = Star_class.star_model(name)
    star.pop_plot_star(0, min_lum_angle, (0.1, 180), f'Figures/Stars/min_lum_{name}.png', pixel_count=200)

