# %%
from Star_class import star_model
import numpy as np 
import matplotlib.pyplot as plt
import os
import astropy.constants as c
import pandas as pd
import matplotlib.colors as mcolors
import matplotlib as mpl
from tqdm import tqdm

mpl.rcParams['font.size'] = 14

# retrieve all the stars in the data directory
# The if statement is there because there is folder created `by MacOS that I cant remove so it doesn't use that.
names = [i for i in os.listdir(os.environ['FMPdata']) if i[0] != '.']
star_data = pd.read_csv('/Users/tobias/Data/Extended_stellar_params.csv', dtype=None, delimiter=';', encoding='utf-8')

# %% 
from pyvo import registry
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
print(np.min(dfw['R*']), np.max(dfw['R*']))
Rmask = (dfw['R*'] > 0.5) * (dfw['R*'] <= 1.0)
df = dfw[Rmask]

# %%
# Creating models and extracting fluxes. This might take a while, so go get a cup of tea and a good night sleep if its on linear interpolation 
def dict_to_array(d):
    new_d = {}
    for key in d.keys():  
        new_d[key] = np.array(d[key])
    return new_d

stars = {'name': [], 'HARDx': [], 'ROSATx': [], 'EUV': [], 'kpr': [], 'Lbol': [], 'association': [], 'Prot': []}
for name in tqdm(names):
    star = star_model(name, interpolation='nearest', verbose=False)
    if "sun" in name.lower(): # Sun data is not found in star_data so we make an exception for those 
        lum_bol = 3.86e33
        association = 'Sun'
        rossby = 2
    else:
        star_idx= star_data[star_data['Star'].str.replace(' ', '') == name.replace('1x-', '')].index
        lum_bol = star_data['Lum'][star_idx].values[0] * 3.86e33 # erg / s
        association = star_data['Associations'][star_idx].values[0]

    x_lum_hard = star.lum_x(image_radius=250, pixel_count=100, wvl_bin=(0.1, 5), grid_type='segmented', nseg=1)
    x_lum_ROSAT = star.lum_x(image_radius=250, pixel_count=100, wvl_bin=(5, 120), grid_type='segmented', nseg=1)
    x_lum_EUV = star.lum_x(image_radius=250, pixel_count=100, wvl_bin=(120, 180), grid_type='segmented', nseg=1)
    kpr_value = np.log10(1.86e-3 * star.params['RotationPeriodStar'] ** -2 * (star.params['RadiusStar'] / c.R_sun.to('cm').value) **-4)

    stars['name'].append(name)
    stars['HARDx'].append(x_lum_hard)
    stars['ROSATx'].append(x_lum_ROSAT)
    stars['EUV'].append(x_lum_EUV)
    stars['kpr'].append(kpr_value)
    stars['Lbol'].append(lum_bol)
    stars['association'].append(association)
    stars['Prot'].append(star.params['RotationPeriodStar'])

stars = dict_to_array(stars)
# %%
# -- Plotting --
associations = np.unique(stars['association'])
fig_pl, ax = plt.subplots(1,1,figsize=(8,8))
fig_kpr, ax2 = plt.subplots(1,3, figsize=(16,5), sharey=True)
fig_lit, ax3 = plt.subplots(1,1, figsize=(8,8))

ax.scatter(np.log10(df['Prot']), df['Lx/bol'], color='k', s=1)
ax3.scatter(np.log10(1.86e-3 * df['R*']**-4 *df['Prot'] **-2), df['Lx/bol'], c='gray', s=2, label='Reiners+2014')

for star in stars['name']:
    idx = np.where(stars['name'] == star)
    ax.scatter(np.log10(stars['Prot'][idx]), np.log10((stars['HARDx'][idx] + stars['ROSATx'][idx] + stars['EUV'][idx])/stars['Lbol'][idx]))

    ax2[2].scatter(stars['kpr'][idx], np.log10(stars['EUV'][idx]/stars['Lbol'][idx]), marker='*', label=stars['association'][idx],
                   c=list(mcolors.TABLEAU_COLORS)[np.where(associations == stars['association'][idx])[0][0]])

    ax2[1].scatter(stars['kpr'][idx], np.log10(stars['ROSATx'][idx]/stars['Lbol'][idx]), marker='*', label=stars['association'][idx],
                   c=list(mcolors.TABLEAU_COLORS)[np.where(associations == stars['association'][idx])[0][0]])
    
    ax2[0].scatter(stars['kpr'][idx], np.log10(stars['HARDx'][idx]/stars['Lbol'][idx]), marker='*', label=stars['association'][idx],
                   c=list(mcolors.TABLEAU_COLORS)[np.where(associations == stars['association'][idx])[0][0]])
    
    ax3.scatter(stars['kpr'][idx], np.log10((10 ** 4) * stars['ROSATx'][idx]/stars['Lbol'][idx]), marker='*', s=100, label=stars['association'][idx][0],
                c=list(mcolors.TABLEAU_COLORS)[np.where(associations == stars['association'][idx])[0][0]])

ax.set_ylabel('Log L$_x$ (erg s$^{-1}$)')
ax.set_xlabel('log P$_{rot}$ (Days)')

ax2[0].invert_xaxis() ; ax2[1].invert_xaxis(); ax2[2].invert_xaxis()
ax2[0].set_ylabel('Log L$_x$ / L$_{bol}$ (erg s$^{-1}$)')
ax2[0].axvline(-3.14, ls='--', color='grey') ; ax2[1].axvline(-3.14, ls='--', color='grey') ; ax2[2].axvline(-3.14, ls='--', color='grey')
ax2[0].set_xlabel('k P$^{-2}$ R$^{-4}$') ; ax2[1].set_xlabel('k P$^{-2}$ R$^{-4}$') ; ax2[2].set_xlabel('k P$^{-2}$ R$^{-4}$')
ax2[0].set_title('Hard Xrays [0.1-5] $\AA$') ; ax2[1].set_title('ROSAT [5-120] $\AA$') ; ax2[2].set_title('EUV [120-180] $\AA$')
fig_kpr.subplots_adjust(wspace=0)

ax3.invert_xaxis()
ax3.set_ylabel('Log L$_x$ / L$_{bol}$ (erg s$^{-1}$)')
ax3.set_xlabel('k P$^{-2}$ R$^{-4}$')
ax3.set_title('ROSAT [5-120] $\AA$')

kpr_axs = fig_lit.gca()
handles, labels = kpr_axs.get_legend_handles_labels()
by_label = dict(zip(labels, handles))
ax3.legend(by_label.values(), by_label.keys())


# fig_pl.savefig('Figures/Lx_Prot_+lit.png', dpi=500, bbox_inches='tight')
# fig_kpr.savefig('Figures/Lx_kpr_Mulitple_bands.pdf', dpi=50,  bbox_inches='tight')
# fig_lit.savefig('Figures/Lx_kpr+Reiners.png', dpi=50, bbox_inches='tight')  

# %%
# ------------------ Fitting ------------------
import pwlf
from scipy.optimize import curve_fit

def pw_linear(x, a, b):
    return np.piecewise(x, [x > -3.14, x <= -3.14], [ lambda x: b, lambda x: a*x + (b - a*(-3.14))])



sat_mask = stars['kpr'] <= -3.14
opt_pw, cov_pw = curve_fit(pw_linear, stars['kpr'], np.log10(stars['ROSATx']/stars['Lbol']))
opt_reiners, cov_reiners = curve_fit(pw_linear, np.log10(1.86e-3 * df['R*']**-4 *df['Prot'] **-2), df['Lx/bol'])
# We expect a unity slope from Reiners
print(opt_pw, np.sqrt(np.diag(cov_pw)))

x_dum = np.linspace(-6, 2, 20)

# fig_fit, ax_fit = plt.subplots(1,1, figsize=(8,8))
# ax_fit.scatter(np.log10(1.86e-3 * df['R*']**-4 *df['Prot'] **-2), df['Lx/bol'], c='gray', s=2, label='Reiners+2014')
ax3.plot(x_dum, 4.5 + pw_linear(x_dum, *opt_pw), c='k', label='Pw fit')
# ax_fit.scatter(stars['kpr'], np.log10(stars['ROSATx']/stars['Lbol']), c='k', label='Data')
# ax_fit.invert_xaxis()
kpr_axs = fig_lit.gca()
handles, labels = kpr_axs.get_legend_handles_labels()
by_label = dict(zip(labels, handles))
ax3.legend(by_label.values(), by_label.keys())

fig_lit.savefig('Figures/Lx_kpr+Reiners+fit.png', dpi=50, bbox_inches='tight')  
plt.show()




# %% 
# --------------- L_x/ AGE plots --------------------

fig_la, ax = plt.subplots(1,1,figsize=(8,8))
sol_lum = 3.86e33

groups = {}
i = 1
for name in names:
    star = star_model(name)
    if "sun" in name.lower(): # Sun data is not found in star_data so we make an exception for those 
        lum_bol = sol_lum
        association = 'Sun'
        rossby = 2
        age_star = 4.5e3
        age_err = 0
    else:
        star_idx= star_data[star_data['Star'].str.replace(' ', '') == name.replace('1x-', '')].index
        lum_bol = star_data['Lum'][star_idx].values[0] * sol_lum # erg / s
        association = star_data['Associations'][star_idx].values[0]
        age_star = star_data['Age'][star_idx].values[0]
        age_err = star_data['Age_err'][star_idx].values[0]
    x_lum = star.lum_x()
    if association not in groups.keys():
        groups[association] = list(mcolors.TABLEAU_COLORS.keys())[i]
        i += 1
        ax.scatter(np.log10(age_star), np.log10(x_lum/lum_bol),
                   color = groups[association], label=association)
    else:
        try:
            ax.scatter(np.log10(age_star), np.log10(x_lum/lum_bol),
                    color = groups[association])
            os.system('clear') # To clear the terminal
        except:
            print(f"Not plotted {name}")

ax.set_ylabel('Log L$_x$ (erg s$^{-1}$)')
ax.set_xlabel('log Age (Myr)')
ax.legend()
fig_la.savefig('Figures/Lx_age_retation.png', dpi=500, bbox_inches='tight')
plt.show()
# %%
# -------------- Lx/open field lines plot -------------------
groups = np.unique(star_data['Associations'].values)

for i, name in enumerate(names):
    if 'sun' in name.lower():
        continue
    star = star_model(name)
    lum_x = star.lum_x(pixel_count=100)
    star_idx = star_data[star_data['Star'].str.replace(' ', '') == name.replace('1x-', '')].index
    frac_open_lines = star_data['phi_open'][star_idx].values[0]
    association = star_data['Associations'][star_idx].values[0]
    p = plt.scatter(frac_open_lines, np.log10(lum_x), label=association, 
                    c=list(mcolors.TABLEAU_COLORS)[np.where(groups == association)[0][0]])


handles, labels = plt.gca().get_legend_handles_labels()
by_label = dict(zip(labels, handles))
plt.xlabel('Fraction of open field lines')
plt.ylabel('log$_{10}$ L$_x$ erg s$^{-1}$')
plt.legend(by_label.values(), by_label.keys())
plt.savefig('Figures/xlum-openfieldlines_wholespec.png', dpi=500, bbox_inches='tight')
plt.show()