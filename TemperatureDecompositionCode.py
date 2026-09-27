import glob
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import salem
import rioxarray
import warnings

plt.rcParams["font.family"] = "Helvetica"
plt.rcParams["axes.titlesize"] = 14
plt.rcParams["axes.labelsize"] = 14
plt.rcParams["xtick.labelsize"] = 12
plt.rcParams["ytick.labelsize"] = 12
plt.rcParams["legend.fontsize"] = 11
plt.rcParams["axes.linewidth"] = 0.8
plt.rcParams["axes.labelweight"] = 'bold'
plt.rcParams["xtick.direction"] = 'out'
plt.rcParams["ytick.direction"] = 'out'
warnings.filterwarnings('ignore')

shp_path = '/Users/zhangperry/Downloads/EU_Shapefile_dissolve.shp'
sigma = 5.670374e-8
bar_height = 0.3


HIST_DIR = '/Volumes/ORCHIDEE/Final_Holiforma_Results/INFORMA.BAU.SSP126.REPEAT/MO/'
BAU_DIR  = '/Volumes/ORCHIDEE/Final_Holiforma_Results/INFORMA.OPTIM.2.1.SSP126.REPEAT/MO/'

hist_srf_files = sorted(glob.glob(HIST_DIR + 'INFORMA.BAU.SSP126.REPEAT_*_1M_sechiba_history.nc'))[-15:]
hist_atm_files = sorted(glob.glob(HIST_DIR + 'INFORMA.BAU.SSP126.REPEAT_*_1M_histmth.nc'))[-15:]
bau_srf_files  = sorted(glob.glob(BAU_DIR  + 'INFORMA.OPTIM.2.1.SSP126.REPEAT_*_1M_sechiba_history.nc'))[-15:]
bau_atm_files  = sorted(glob.glob(BAU_DIR  + 'INFORMA.OPTIM.2.1.SSP126.REPEAT_*_1M_histmth.nc'))[-15:]

# for label, files in [('HIST SRF', hist_srf_files), ('HIST ATM', hist_atm_files),
#                       ('BAU SRF',  bau_srf_files),  ('BAU ATM',  bau_atm_files)]:
#     for f in files: print("  ", f)
    
def load_multiyear_srf(files):
    if not files:
        raise FileNotFoundError(f"No files. Please check.")
    
    ds = xr.open_mfdataset(
        files,
        combine='nested',
        concat_dim='time_counter',
        engine='netcdf4'
    )
    ds = ds.isel(lat=slice(13, 83), lon=slice(50, 94))
    ds['G']     = ds['Qg']
    ds['tsol']  = ds['temp_sol'] + 273.15
    ds['eps_a'] = ds['lwdown'] / (sigma * (ds['tair']**4))
    subset = ds.mean(dim='time_counter')
    return subset

def load_multiyear_atm(files, swdown):
    if not files:
        raise FileNotFoundError(f"No files. Please check.")
    ds = xr.open_mfdataset(
        files,
        combine='nested',
        concat_dim='time_counter',
        engine='netcdf4')
    ds = ds.isel(lat=slice(13, 83), lon=slice(50, 94))
    ds['alb_avg'] = ds['SWupSFC'] / ds['SWdnSFC']
    subset = ds.mean(dim='time_counter')
    return subset


ds_hist_srf = load_multiyear_srf(hist_srf_files)
ds_bau_srf  = load_multiyear_srf(bau_srf_files)

ds_hist = load_multiyear_atm(hist_atm_files, ds_hist_srf['swdown'])
ds_bau  = load_multiyear_atm(bau_atm_files,  ds_bau_srf['swdown'])


delta = ds_bau - ds_hist
delta_srf = ds_bau_srf - ds_hist_srf

beta = 1.0 / (4 * sigma * ds_hist_srf['eps_a'] * ds_hist_srf['tair']**3)

contributions = {
    'alb':   beta * (ds_hist_srf['swdown'] * delta['alb_avg']),
    'Rsi':   beta * (-(1 - ds_hist['alb_avg']) * delta_srf['swdown']),
    'le':    beta * delta_srf['fluxlat'],
    'h':     beta * delta_srf['fluxsens'],
    'g':     beta * delta_srf['G'],
    'eps_a': beta * (- sigma * ds_hist_srf['tair']**4 * delta_srf['eps_a']),
    'ts':    beta * (4 * sigma * (ds_hist_srf['tsol']**3) * delta_srf['tsol'])}

residual = sum(contributions.values()) - delta_srf['tair']
shdf = salem.read_shapefile(shp_path)
def get_clipped_zonal(da):
    return da.rio.write_crs("EPSG:4326").rio.clip(shdf.geometry.values, shdf.crs, all_touched=True).mean(dim='lon')
zonal_data = {
    '$\\Delta T_a|\\epsilon$': get_clipped_zonal(contributions['eps_a']),
    '$\\Delta T_a|G$':         get_clipped_zonal(contributions['g']),
    '$\\Delta T_a|LE+H$':      get_clipped_zonal(contributions['le'] + contributions['h']),
    '$\\Delta T_a|R_{si}$':    get_clipped_zonal(contributions['Rsi']),
    '$\\Delta T_a|\\alpha$':   get_clipped_zonal(contributions['alb']),
    '$\\Delta T_a|T_s$':       get_clipped_zonal(contributions['ts'])}

residual_zonal = get_clipped_zonal(residual)
total_ta_zonal = get_clipped_zonal(delta_srf['tair'])
ta_decomposition_zonal = sum(contributions.values()).rio.write_crs("EPSG:4326").rio.clip(shdf.geometry.values, shdf.crs, all_touched=True).mean(dim='lon')


lats = total_ta_zonal.lat.values
colors = ['#9e2a2b', '#fc8d62', '#8da0cb', '#e78ac3', '#a6d854', '#ffd92f']


fig, axes = plt.subplots(2, 2, figsize=(14, 12), gridspec_kw={'hspace': 0.2, 'wspace': 0.1, 'width_ratios': [2, 1]})
def clip_and_plot_map(da, ax, title, cmap, vmin, vmax):
    da_clipped = da.rio.write_crs("EPSG:4326").rio.clip(shdf.geometry.values, shdf.crs, all_touched=True)
    im = da_clipped.plot(ax=ax, cmap=cmap, vmin=vmin, vmax=vmax, add_colorbar=True,
                         cbar_kwargs={'label': 'K', 'shrink': 0.8})
    shdf.boundary.plot(ax=ax, color='black', linewidth=0.8)
    ax.set_title(title, fontweight='bold')
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")

clip_and_plot_map(delta_srf['tair'], axes[0, 0], r"$\Delta T_a^{\!\mathrm{Actual}}$ (Annual)", 'RdBu_r', -2, 2)
clip_and_plot_map(residual, axes[1, 0], r"$\Delta T_a^{\!\mathrm{Actual}}$ - $\Delta T_a^{\!\mathrm{Decomposed}}$ (Annual)", 'PRGn', -4, 4)
pos_base = np.zeros(len(lats))
neg_base = np.zeros(len(lats))

for i, (name, zonal_vals) in enumerate(zonal_data.items()):
    vals = np.nan_to_num(zonal_vals.values)
    starts = np.where(vals >= 0, pos_base, neg_base)
    axes[0, 1].barh(lats, vals, left=starts, height=bar_height, 
                    label=name, color=colors[i], alpha=0.8, edgecolor='none')
    pos_base += np.maximum(vals, 0)
    neg_base += np.minimum(vals, 0)

axes[0, 1].plot(ta_decomposition_zonal, lats, color='red', lw=1.5, label=r'$\Delta T_a^{\!\mathrm{Decomposed}}$')
axes[0, 1].plot(total_ta_zonal, lats, color='black', lw=1.5, label=r'$\Delta T_a^{\!\mathrm{Actual}}$')
axes[0, 1].axvline(0, color='black', lw=0.8, linestyle='-')
axes[0, 1].set_xlabel("$\Delta T_a$ (K)")
axes[0, 1].set_ylabel("Latitude")
axes[0, 1].set_xlim(-4, 4)
axes[0, 1].set_ylim(35, 71)
axes[0, 1].grid(axis='x', linestyle='--', alpha=0.3)
axes[0, 1].legend(loc='lower center', bbox_to_anchor=(0.23, 0.0), ncol=1, frameon=False)

axes[1, 1].plot(residual_zonal, lats, color='#5aae61', lw=2, label='Zonal Residual')
axes[1, 1].axvline(0, color='black', lw=0.8, linestyle='-')
axes[1, 1].set_xlabel("Residual (K)")
axes[1, 1].set_ylabel("Latitude")
axes[1, 1].set_xlim(-4, 4)
axes[1, 1].set_xticks([-4, -2, 0, 2, 4])
axes[1, 1].set_xticklabels(['-4', '-2', '0', '2', '4'])
axes[1, 1].set_ylim(35, 71)
axes[1, 1].grid(axis='both', linestyle='--', alpha=0.3)
axes[1, 1].fill_betweenx(lats, 0, residual_zonal, color='#5aae61', alpha=0.2)
plt.show()