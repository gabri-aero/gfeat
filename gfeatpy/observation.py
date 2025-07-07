from ._core.observation import *


from gfeatpy._core import observation as _observation

__all__ = []

for name in dir(_observation):
    obj = getattr(_observation, name)
    if getattr(obj, "__module__", "").startswith("gfeatpy._core.observation"):
        try:
            obj.__module__ = "gfeatpy.observation"
        except AttributeError:
            pass  # some objects like modules or builtins might not support this
        globals()[name] = obj
        __all__.append(name)