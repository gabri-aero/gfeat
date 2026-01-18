import importlib

_core = importlib.import_module("gfeatpy._core")

globals().update(
    {name: getattr(_core, name) for name in dir(_core) if not name.startswith("_")}
)

del _core
del importlib

# Define what should be visible to users
__all__ = [name for name in globals() if not name.startswith("_")]