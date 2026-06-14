import os
import shutil
import subprocess
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

PACKAGES = [
    "maass_forms_klein",
    "maass_forms_klein.functions",
    "maass_forms_klein.database",
    "maass_forms_klein.utils",
    "maass_forms_klein.hyperbolic_space",
    "maass_forms_klein.modform",
    "maass_forms_klein.search",
]

if HAS_SAGE:
    from Cython.Build import cythonize
    from setuptools.extension import Extension

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

    extra_compile_args = [
        "-Wno-unused-function",
        "-Wno-implicit-function-declaration",
        "-Wno-unused-variable",
        "-Wno-deprecated-declarations",
        "-Wno-deprecated-register",
    ]

    ext_modules = [
        Extension(
            "maass_forms_klein.hyperbolic_space.upper_half_space",
            sources=[os.path.join("src/maass_forms_klein/hyperbolic_space/upper_half_space.pyx")],
            extra_compile_args=extra_compile_args,
            include_dirs=INCLUDE_DIRS,
            library_dirs=LIBRARY_DIRS,
        )
    ]

    try:
        from sage.misc.package_dir import cython_namespace_package_support
    except ImportError:

        @contextmanager
        def cython_namespace_package_support():
            yield

    with cython_namespace_package_support():
        setuptools.setup(
            packages=PACKAGES,
            ext_modules=cythonize(
                ext_modules,
                include_path=["src", SAGE_LIB],
                compiler_directives={
                    "embedsignature": True,
                    "language_level": "3",
                },
            ),
        )
else:
    import warnings

    warnings.warn(
        "SageMath/PassageMath not found. Installing without Cython extensions.",
        stacklevel=1,
    )
    setuptools.setup(packages=PACKAGES)
