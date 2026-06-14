import os
import shutil
import subprocess
import sys
import sysconfig
from contextlib import contextmanager

import setuptools

# Detect SageMath/PassageMath. Cython extensions require it, but we allow
# installation without it so that the pure-Python parts remain usable.
SAGE_LIB = None
try:
    from sage.env import SAGE_LIB
except ModuleNotFoundError:
    SAGE_LIB = os.getenv("SAGE_LIB") or sysconfig.get_path("purelib")

HAS_SAGE = SAGE_LIB is not None and os.path.isdir(SAGE_LIB + "/sage")

if HAS_SAGE:
    from Cython.Build import cythonize
    from setuptools.extension import Extension

    INCLUDE_DIRS = []
    LIBRARY_DIRS = []

    def _add_prefix(prefix):
        inc = os.path.join(prefix, "include")
        lib = os.path.join(prefix, "lib")
        if os.path.isdir(inc):
            INCLUDE_DIRS.append(inc)
        if os.path.isdir(lib):
            LIBRARY_DIRS.append(lib)

    if shutil.which("brew") is not None:
        try:
            _add_prefix(
                subprocess.check_output(
                    ["brew", "--prefix"], text=True, stderr=subprocess.DEVNULL
                ).strip()
            )
        except (subprocess.CalledProcessError, FileNotFoundError):
            pass

    if os.environ.get("CONDA_PREFIX"):
        _add_prefix(os.environ["CONDA_PREFIX"])

    # pkg-config picks up custom installs (e.g. /opt). Distro packages
    # like libmpc-dev land in /usr/include and are found by gcc directly.
    if shutil.which("pkg-config") is not None:
        for _lib in ("mpc", "mpfr", "gmp", "flint"):
            for _flag, _target in (
                ("--cflags-only-I", INCLUDE_DIRS),
                ("--libs-only-L", LIBRARY_DIRS),
            ):
                try:
                    _out = subprocess.check_output(
                        ["pkg-config", _flag, _lib],
                        text=True,
                        stderr=subprocess.DEVNULL,
                    ).split()
                except subprocess.CalledProcessError:
                    continue
                _target.extend(a[2:] for a in _out if a.startswith(("-I", "-L")))

    INCLUDE_DIRS = list(dict.fromkeys(INCLUDE_DIRS))
    LIBRARY_DIRS = list(dict.fromkeys(LIBRARY_DIRS))

    extra_compile_args = [
        "-Wno-unused-function",
        "-Wno-implicit-function-declaration",
        "-Wno-unused-variable",
        "-Wno-deprecated-declarations",
        "-Wno-deprecated-register",
    ]

    ext_modules = [
        Extension(
            "maass_form_core.functions.bessel.besselk_dp",
            ["src/maass_form_core/functions/bessel/besselk_dp.pyx"],
            include_dirs=INCLUDE_DIRS,
            extra_compile_args=extra_compile_args,
            library_dirs=LIBRARY_DIRS,
        ),
        Extension(
            "maass_form_core.functions.functions_cy",
            ["src/maass_form_core/functions/functions_cy.pyx"],
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
                include_path=["src", SAGE_LIB],
                compiler_directives={
                    "embedsignature": True,
                    "language_level": "3",
                },
                annotate=True,
            ),
        )
else:
    print(
        "WARNING: SageMath/PassageMath not found. Installing without Cython extensions.",
        file=sys.stderr,
    )
    setuptools.setup()
