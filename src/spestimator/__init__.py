try:
    from importlib.metadata import version, PackageNotFoundError

    __version__ = version("spestimator")
except PackageNotFoundError:
    __version__ = "unknown"
