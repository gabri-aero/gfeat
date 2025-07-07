import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm, Normalize
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from gfeatpy.observation import *
from gfeatpy.utils import *


def plot_rms_per_coefficient_per_degree(ax, sh_coefs, l_min=2, show_kaula=False):
    sigma_l = np.sqrt(np.array(sh_coefs.degree_variance()))
    sigma_l[sigma_l <= 0] =np.nan
    l = np.arange(l_min, len(sigma_l))
    sigma_l = sigma_l[l_min:] / np.sqrt(2*l+1)
    ax.plot(l, sigma_l)
    ax.set_yscale('log')
    ax.set_xlabel('Degree')
    ax.set_ylabel('$\\delta_l$')
    if show_kaula:
        ax.plot(l, np.sqrt(160e-12/l**3), 'r--')


def plot_pyramid_coefs(ax, sh_coefs, zonal_terms=False):
    coefs = sh_coefs.get_sigma_x()
    if not zonal_terms:
        coefs[:, 0] = 0
    l_max = coefs.shape[0]-1
    Clm = np.tril(coefs)
    Slm = np.fliplr(np.triu(coefs, k=1).T)
    pyramid = np.hstack((Slm, Clm))
    pyramid[pyramid == 0] = np.nan
    im = ax.imshow(pyramid, cmap='jet', extent=[-l_max, l_max, 0, l_max], origin='lower', norm=LogNorm())
    fig = plt.gcf()
    fig.colorbar(im, ax=ax)
    ax.invert_yaxis()


def plot_summary(sh_coefs, axs):
    plot_rms_per_coefficient_per_degree(axs[0], sh_coefs)
    plot_pyramid_coefs(axs[1], sh_coefs)


def get_basemap_axes():
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.add_feature(cfeature.COASTLINE, linewidth=0.7)
    ax.set_xlabel('Longitude [deg]')
    ax.set_ylabel('Latitude [deg]')
    ax.set_yticks(np.linspace(-90, 90, 7))
    ax.set_xticks(np.linspace(-180, 180, 13))
    return ax


def adjust_colorbar_ticks(cb):
    vmin, vmax = cb.mappable.get_clim()
    rango = vmax - vmin

    # Elegir orden de magnitud para el paso
    orden = 10 ** np.floor(np.log10(rango/4))
    
    # Elegir paso como múltiplo 'perfecto' (1, 2, 5, 10, ...)
    for factor in [1, 2, 5, 10]:
        paso = orden * factor
        num_ticks = (vmax - vmin) / paso
        if 3 <= num_ticks <= 10:
            break

    # Generar ticks dentro del rango
    tick_start = np.ceil(vmin / paso) * paso
    tick_end   = np.floor(vmax / paso) * paso
    ticks = np.arange(tick_start, tick_end + paso/2, paso)

    # Aplicar los ticks
    cb.set_ticks(ticks)
    cb.set_ticklabels([f'{t:g}' for t in ticks])


def base_synthesis(lon, lat, z, functional, ax=None):
    # Define plot settings based on synthesis functional type
    if ax is None: 
        fig = plt.figure(figsize=(7,3.5))
        ax = get_basemap_axes()

    zlim = np.max(np.abs(z))
        
    if type(functional).__name__ == "GravityAnomaly":
        name = "$\Delta$g"
        unit = "[mGal]"
    elif type(functional).__name__ == "EquivalentWaterHeight":
        name = "EWH"
        if zlim < 0.2:
            unit = "[mm]"
            z = z*1e3
            zlim = zlim*1e3
        elif zlim < 2:
            unit = "[cm]"
            z = z*1e2
            zlim = zlim*1e2
        else:
            unit = "[m]"
    elif type(functional).__name__ == "GeoidHeight":
        name = "N"
        unit = "[m]"

    if np.min(z) < 0:
        levels = np.linspace(-zlim, zlim, 101)
        clabel = f"{name} {unit}"
    else:
        levels = np.linspace(0, zlim, 101)
        clabel = f"$\sigma$({name}) {unit}"

    contour = ax.contourf(lon, lat, z, levels=levels, alpha=0.6, cmap='seismic', antialiased=True)
    cbar = plt.colorbar(contour, ax=ax, shrink=0.5)
    cbar.set_label(clabel)
    adjust_colorbar_ticks(cbar)
    plt.tight_layout()
    return contour, cbar


def synthesis(sh_coefs, n_lon, n_lat, functional,  ax=None):
    [lon, lat, z] = sh_coefs.synthesis(n_lon, n_lat, functional)
    return base_synthesis(lon, lat, z, functional, ax)


def plot_ground_track(I, Nr, Nd, we_0, ax, *args, **kwargs):
    Tr = 1
    wo_dot = 2*np.pi * Nr / Tr
    we_dot = -2*np.pi * Nd / Tr
    t = np.linspace(0, Tr, 100_000)
    
    # Definition of wo, we
    wo_0 = 0  # for simplicity
    wo = wo_0 + wo_dot * t
    we = we_0 + we_dot * t
    
    # Conversion to x, y, z
    x = np.cos(we) * np.cos(wo) - np.sin(we) * np.sin(wo) * np.cos(I)
    y = np.sin(we) * np.cos(wo) + np.cos(we) * np.sin(wo) * np.cos(I)
    z = np.sin(wo) * np.sin(I)
    
    # Conversion to lon, lat
    lon = np.arctan2(y, x)
    lat = np.arcsin(z)

    # Convert
    lon = np.rad2deg(lon)
    lat = np.rad2deg(lat)
    
    # Plot ground track
    ax.scatter(lon, lat, *args, **kwargs)
    
    
show = plt.show
savefig = plt.savefig
gca = plt.gca
figure = plt.figure