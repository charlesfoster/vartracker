"""Version information for the vartracker package."""

from __future__ import annotations

try:  # Python 3.8+
    from importlib import metadata
except ImportError:  # pragma: no cover -- fallback for very old Pythons
    import importlib_metadata as metadata  # type: ignore

# NOTE: When bumping the project version remember to update this fallback value
# alongside the version declared in pyproject.toml.
_FALLBACK_VERSION = "1.0.1"

try:
    __version__ = metadata.version("vartracker")
except metadata.PackageNotFoundError:  # pragma: no cover - during local dev
    __version__ = _FALLBACK_VERSION

__all__ = ["__version__"]
