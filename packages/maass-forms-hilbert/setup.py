import os
import shutil
import subprocess
import sysconfig
from contextlib import contextmanager

import setuptools
from Cython.Build import cythonize
from setuptools.extension import Extension

# Find correct value for SAGE_LIB which is needed to compile the Cython extensions.
try:
    from sage.env import SAGE_LIB
except ModuleNotFoundError:
    SAGE_LIB = os.getenv("SAGE_LIB") or sysconfig.get_path("purelib")
# Check if sage is installed in SAGE_LIB
if not os.path.isdir(SAGE_LIB + "/sage"):
    raise ModuleNotFoundError("Sagemath or passagemath needs to be installed.")

INCLUDE_DIRS = []
LIBRARY_DIRS = []
if shutil.which("brew") is not None:
    proc = subprocess.Popen(
        "brew --prefix",
        shell=True,
        stdout=subprocess.PIPE,
        stdin=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        close_fds=True,
    )
    HOMEBREW_PREFIX = proc.stdout.readline().decode("utf-8").strip()
    HOMEBREW_LIB = HOMEBREW_PREFIX + "/lib"
    LIBRARY_DIRS.append(HOMEBREW_LIB)
    HOMEBREW_INC = HOMEBREW_PREFIX + "/include"
    INCLUDE_DIRS.append(HOMEBREW_INC)

# Extension modules using Cython
extra_compile_args = [
    "-Wno-unused-function",
    "-Wno-implicit-function-declaration",
    "-Wno-unused-variable",
    "-Wno-deprecated-declarations",
    "-Wno-deprecated-register",
]

# Find maass_form_core source for cimport resolution
try:
    import maass_form_core
    CORE_SRC = os.path.dirname(os.path.dirname(maass_form_core.__file__))
except ImportError:
    CORE_SRC = os.path.join(os.path.dirname(__file__), "..", "maass-form-core", "src")

ext_modules = [
    Extension(
        "maass_forms_hilbert.functions.functions_cy",
        ["src/maass_forms_hilbert/functions/functions_cy.pyx"],
        include_dirs=INCLUDE_DIRS,
        extra_compile_args=extra_compile_args,
        library_dirs=LIBRARY_DIRS,
    ),
]

try:
    from sage.misc.package_dir import cython_namespace_package_support
except ImportError:

    @contextmanager
    def cython_namespace_package_support():
        yield


with cython_namespace_package_support():
    setuptools.setup(
        ext_modules=cythonize(
            ext_modules,
            include_path=["src", SAGE_LIB, CORE_SRC],
            compiler_directives={
                "embedsignature": True,
                "language_level": "3",
            },
            annotate=True,
        ),
    )
