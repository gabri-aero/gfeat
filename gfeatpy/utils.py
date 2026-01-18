from ._core.utils import *

from gfeatpy._core import utils as _utils

__all__ = []

for name in dir(_utils):
    obj = getattr(_utils, name)
    if getattr(obj, "__module__", "").startswith("gfeatpy._core.utils"):
        try:
            obj.__module__ = "gfeatpy.utils"
        except AttributeError:
            pass  # some objects like modules or builtins might not support this
        globals()[name] = obj
        __all__.append(name)