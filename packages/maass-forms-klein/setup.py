import os
import shutil
import subprocess
import sys

import setuptools
from sage_setup.extensions import create_extension
from setuptools.extension import Extension
import Cython.Compiler.Main
from Cython.Build import cythonize
from sage.env import SAGE_LIB

debug = False
gdb_debug = True
if os.environ.get('SAGE_DEBUG', None) == 'yes':
    print('Enabling Cython debugging support')
    debug = True
    Cython.Compiler.Main.default_options['gdb_debug'] = True
    Cython.Compiler.Main.default_options['output_dir'] = 'build'
    gdb_debug = True

LIBRARY_DIRS = []
INCLUDE_DIRS = []
if shutil.which('brew') is not None:
    proc = subprocess.Popen("/opt/homebrew/bin/brew --prefix", shell=True,
                            stdout=subprocess.PIPE, stdin=subprocess.PIPE,
                            stderr=subprocess.STDOUT, close_fds=True)
    HOMEBREW_PREFIX = proc.stdout.readline().decode('utf-8').strip()
    HOMEBREW_LIB = HOMEBREW_PREFIX + '/lib'
    LIBRARY_DIRS.append(HOMEBREW_LIB)
    HOMEBREW_INC = HOMEBREW_PREFIX + '/include'
    INCLUDE_DIRS.append(HOMEBREW_INC)

INCLUDE_DIRS += ['src']
extra_compile_args = ['-Wno-unused-function',
                      '-Wno-implicit-function-declaration',
                      '-Wno-unused-variable',
                      '-Wno-deprecated-declarations',
                      '-Wno-deprecated-register',
                      '-Wno-unreachable-code',
                      '-Wno-unreachable-code-fallthrough']

ext_modules = [
    Extension(
        'knot_maass.functions.besselk_dp',
        sources=[os.path.join('src/knot_maass/functions/besselk_dp.pyx')],
        extra_compile_args=extra_compile_args,
        include_dirs=INCLUDE_DIRS, library_dirs=LIBRARY_DIRS
    ),
    Extension(
        'knot_maass.hyperbolic_space.upper_half_space',
        sources=[os.path.join('src/knot_maass/hyperbolic_space/upper_half_space.pyx')],
        extra_compile_args=extra_compile_args,
        include_dirs=INCLUDE_DIRS, library_dirs=LIBRARY_DIRS
    )
]

print("ext modules=",ext_modules)
setuptools.setup(
    packages=['knot_maass',
              'knot_maass.functions',
              'knot_maass.hyperbolic_space',
              'knot_maass.modform'],
    ext_modules=cythonize(
        ext_modules,
        include_path=['src', SAGE_LIB],
        compiler_directives={
            'embedsignature': True,
            'language_level': '3',
        },
    ),
)
