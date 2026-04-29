import logging
from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("maass_forms_hilbert")
except PackageNotFoundError:
    __version__ = "0.0.0"

logging.basicConfig(level=logging.ERROR, format="%(message)s")
logging.captureWarnings(True)
log = logging.getLogger(__name__)
