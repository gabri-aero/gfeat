import matplotlib.pyplot as plt
from matplotlib.axes import Axes
import numpy as np
from matplotlib.colors import LogNorm, Normalize
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from gfeatpy.gravity import BaseFunctional
from gfeatpy.utils import ma_dot, raan_dot, aop_dot, sma_from_rgt
from gfeatpy import planet
    

def pyramid(base_object, colormap='jet'):
    """
    Function to plot the SH coefficients in a pyramid format.
    
    :param base_object: Any object containing SH coefficients or errors
    :param colormap: Colormap to use for plotting
    """
    coefs = base_object.get_sigma_x()
    l_max = coefs.shape[0]-1
    Clm = np.tril(coefs)
    Slm = np.fliplr(np.triu(coefs, k=1).T)
    pyramid = np.hstack((Slm, Clm))
    im = plt.imshow(pyramid, cmap=colormap, extent=[-l_max, l_max, 0, l_max], origin='lower', norm=LogNorm())
    ax = plt.gca()
    ax.invert_yaxis()
    plt.colorbar(im, ax=ax)
    return im
    
    
def synthesis(base_object, n_lon: int, n_lat: int, functional: BaseFunctional, colormap = 'seismic', z_max: float = None):
    """
    Plot the global SH synthesis on an Earth map.
    (To DO: add option to remove Earth map background)
    
    :param base_object: Any object containing SH coefficients or errors
    :param n_lon: Number of longitude points
    :type n_lon: int
    :param n_lat: Number of latitude points
    :type n_lat: int
    :param functional: Any synthesis functional: GravityAnomaly, EquivalentWaterHeight, GeoidHeight
    :type functional: BaseFunctional
    :param colormap: Colormap to use for plotting
    :type colormap: str
    :param z_max: Cutting value for the color scale
    :type z_max: float
    """
    def adjust_colorbar_ticks(cb):
        # Define range
        vmin, vmax = cb.mappable.get_clim()
        range = vmax - vmin
        # Set order of magnitudes
        order = 10 ** np.floor(np.log10(range/4))
        # Set step as an integer
        for factor in [1, 2, 5, 10]:
            step = order * factor
            num_ticks = (vmax - vmin) / step
            if 3 <= num_ticks <= 10:
                break
        # Generate ticks
        tick_start = np.ceil(vmin / step) * step
        tick_end   = np.floor(vmax / step) * step
        ticks = np.arange(tick_start, tick_end + step/2, step)
        # Apply ticks
        cb.set_ticks(ticks)
        ticks[np.abs(ticks) < 1e-10] = 0
        cb.set_ticklabels([f'{t:g}' for t in ticks])


    def base_synthesis(lon, lat, z, functional, z_max):
        # Define plot settings based on synthesis functional type
        if z_max is None:
            z_max = np.max(np.abs(z))
        # Clip values to +/- z_max
        z[z > z_max] = z_max
        z[z < -z_max] = -z_max 
        # Automatically define unit based on functional, make it also pretty
        if type(functional).__name__ == "GravityAnomaly":
            name = "$\\Delta$g"
            unit = "[mGal]"
        elif type(functional).__name__ == "EquivalentWaterHeight":
            name = "EWH"
            if z_max < 0.2:
                unit = "[mm]"
                z = z*1e3
                z_max = z_max*1e3
            elif z_max < 2:
                unit = "[cm]"
                z = z*1e2
                z_max = z_max*1e2
            else:
                unit = "[m]"
        elif type(functional).__name__ == "GeoidHeight":
            name = "N"
            unit = "[m]"
        # Define levels for contour plot and colorbar units
        if np.min(z) < 0:
            levels = np.linspace(-z_max, z_max, 101)
            clabel = f"{name} {unit}"
        else:
            levels = np.linspace(0, z_max, 101)
            clabel = f"$\\sigma$({name}) {unit}"
        # Plot contour field
        ax = plt.gca()
        contour = ax.contourf(lon, lat, z, levels=levels, alpha=0.6, cmap=colormap, antialiased=True)
        # Set up colorbar
        cbar = plt.colorbar(contour, ax=ax, shrink=0.5)
        cbar.set_label(clabel)
        adjust_colorbar_ticks(cbar)
        # For nice display -> tight layout
        plt.tight_layout()
        # Return contour and colorbar objects to allow for further modification
        return contour, cbar
    
    def get_basemap_axes():
        # Define projection
        ax = plt.axes(projection=ccrs.PlateCarree())
        # Add continents coastline as background
        ax.add_feature(cfeature.COASTLINE, linewidth=0.7)
        # Set up axes
        ax.set_xlabel('Longitude [deg]')
        ax.set_ylabel('Latitude [deg]')
        ax.set_yticks(np.linspace(-90, 90, 7))
        ax.set_xticks(np.linspace(-180, 180, 13))
        return ax

    get_basemap_axes()
    [lon, lat, z] = base_object.synthesis(n_lon, n_lat, functional)
    return base_synthesis(lon, lat, z, functional, z_max)


def ground_track(Nr, Nd, I, we_0 = 0, wo_0 = 0, samples_per_rev = 3000, **kwargs):
    """
    Function to plot any repeating ground track. Generally to be used along with :func:`gfeatpy.plotting.synthesis` to observe ground-track correlation

    :param Nr: Number of revolutions
    :param Nd: Nodal days
    :param I: Inclination [rad]
    :param we_0: Initial Earth-fixed longitude of ascending node [rad]
    :param wo_0: Initial Earth-fixed argument of latitude [rad]
    :param samples_per_rev: Samples per revolution (higher = smoother ground-track)
    :param kwargs: Additional keyword arguments to be passed to plt.scatter
    """
    Tr = 1
    wo_dot = 2*np.pi * Nr / Tr
    we_dot = -2*np.pi * Nd / Tr
    n_samples = Nr * samples_per_rev
    t = np.linspace(0, Tr, n_samples)
    
    # Definition of wo, we
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
    plt.scatter(lon, lat, **kwargs)


class DegreeAmplitudePlotter:
    """
    Class to plot the RMS per coefficient per degree spectrum

    It allows for adding as many objects as needed with a lot
    of flexibility for formatting
    """
    def __init__(self, figsize, functional: BaseFunctional = None):
        """
        Constructor for DegreeAmplitudePlotter class.

        :param figsize: matplotlib figure size
        :type figsize: tuple
        :param functional: BaseFunctional object to scale be applied to the RMS values
        :type functional: BaseFunctional
        """
        self.f = (
            (lambda l: np.ones_like(l))
            if functional is None 
            else (lambda l: functional._common_factor * np.array([functional._degree_dependent_factor(li) for li in l]))
        )
        self.fig = plt.figure(figsize=figsize)
        self.ax = plt.gca()
        # Set Y axis unit and name
        if type(functional).__name__ == "GravityAnomaly":
            name = "$\\Delta$g"
            unit = "[mGal]"
        elif type(functional).__name__ == "EquivalentWaterHeight":
            name = "EWH"
            unit = "[m]"
        elif type(functional).__name__ == "GeoidHeight":
            name = "N"
            unit = "[m]"
        else:
            name = "$\\delta_l$"
            unit = ""
        self.ax.set_ylabel(f"{name} {unit}")

    def add_item(self, object, *args, show_error = None, **kwargs):
        """
        Add an item to the degree amplitude plot.
        
        :param object: Any object containing SH coefficients or errors
        :param args: Additional positional arguments to be passed to ax.semilogy
        :param show_error: A boolean to indicate whether to plot the error (True) or the signal (False). If None, the signal is plotted by default.
        :param kwargs: Additional keyword arguments to be passed to ax.semilogy
        """
        # Distinguish between plotting error and signal power
        if show_error is None:
            s = object.rms_per_coefficient_per_degree()
        else:
            s = object.rms_per_coefficient_per_degree(show_error)
        l = np.arange(0, len(s))
        s = self.f(l) * s
        s[s < 1e-200] = np.nan
        self.ax.semilogy(l, s, *args, **kwargs)

    def show():
        """
        Wrapper for show method
        """
        plt.show()



class GroundTrackExplorer:
    """
    Class to explore repeating ground tracks for the global planet object as defined in gfeatpy.planet.
    It is essential for orbit selection since Nr and Nd are input arguments for analytical observations
    within gfeatpy.observation module.
    """
    def __init__(self, h_min: float, h_max: float, Nd_max: int):
        """
        Constructor for GroundTrackExplorer class.

        :param h_min: Minimum orbit height [m]
        :param h_max: Maximum orbit height [m]
        :param Nd_max: Maximum number of nodal days to consider
        """
        self.r_min = planet.ae + h_min
        self.r_max = planet.ae + h_max
        self.r_mean = (self.r_min + self.r_max) / 2
        self.Nd_max = Nd_max

    def run(self, I):
        """
        Explore repeating ground tracks within explorer bounds for given inclination.

        :param I: Inclination [rad]
        """
        def f(sma, I):
            we_dot = raan_dot(sma, 0, I) - planet.theta_dot
            wo_dot = aop_dot(sma, 0, I) + ma_dot(sma, 0, I)
            return - wo_dot / we_dot
        # Compute boundaries for fraction Nr/Nd
        day = 2 * np.pi / planet.theta_dot  # day in secs
        fmax = f(self.r_min, I)
        fmin = f(self.r_max, I)
        # Loop through all possibles nodal days
        for nd in np.arange(1, self.Nd_max+1):
            nd = int(nd)
            # Loop through all possible number of revolutions that are within the
            # semi-major axis constraints
            for nr in np.arange(np.ceil(nd*fmin), np.floor(nd*fmax+1)):
                nr = int(nr)
                # Make sure the fraction is not reducible
                if nr % nd != 0 or nd == 1:
                    sma = sma_from_rgt(nr, nd, self.r_mean, I, 0, 100, 1.0, 1e-3)
                    wo_dot = aop_dot(sma, 0, I) + ma_dot(sma, 0, I)
                    t = 2*np.pi*nr / wo_dot / day
                    h = (sma - planet.ae) / 1e3
                    plt.plot(t, h, 'ko')
                    plt.text(t, h, f'{nr}', fontsize=10)
        plt.ylabel('Orbit height [km]')
        plt.xlabel('Repeatability periods [days]')
        plt.ylim([(self.r_min - planet.ae) / 1e3, (self.r_max - planet.ae) / 1e3])

    def show(self):
        """
        Wrapper for show method
        """
        plt.show()


        
    

    