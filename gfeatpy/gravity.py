from ._core.gravity import *
from gfeatpy._core import gravity as _gravity

__all__ = []

for name in dir(_gravity):
    obj = getattr(_gravity, name)
    if getattr(obj, "__module__", "").startswith("gfeatpy._core.gravity"):
        try:
            obj.__module__ = "gfeatpy.gravity"
        except AttributeError:
            pass  # some objects like modules or builtins might not support this
        globals()[name] = obj
        __all__.append(name)


import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm


def pyramid_plot(self, ax, zonal_terms = False):
    coefficients = np.abs(self.coefficients)
    if not zonal_terms:
        coefficients[:, 0] = 0
    l_max = coefficients.shape[0]-1
    Clm = np.tril(coefficients)
    Slm = np.fliplr(np.triu(coefficients, k=1).T)
    pyramid = np.hstack((Slm, Clm))
    pyramid[pyramid == 0] = np.nan
    im = ax.imshow(pyramid, cmap='jet', extent=[-l_max, l_max, 0, l_max], origin='lower', norm=LogNorm())
    fig = plt.gcf()
    fig.colorbar(im, ax=ax)
    ax.invert_yaxis()


SphericalHarmonics.pyramid_plot = pyramid_plot