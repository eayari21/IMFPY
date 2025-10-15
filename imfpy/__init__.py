"""IMFPY simulation and visualization toolkit."""

from importlib.metadata import PackageNotFoundError, version

__all__ = ["__version__"]

try:
    __version__ = version("imfpy")
except PackageNotFoundError:  # pragma: no cover - local source tree
    __version__ = "0.1.dev0"
