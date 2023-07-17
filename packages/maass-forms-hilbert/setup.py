import os
import shutil
import subprocess
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

# INCLUDE_DIRS += ['src/hilbert_maass/modules']
extra_compile_args = ['-Wno-unused-function',
                      '-Wno-implicit-function-declaration',
                      '-Wno-unused-variable',
                      '-Wno-deprecated-declarations',
                      '-Wno-deprecated-register',
                      '-Wno-unreachable-code',
                      '-Wno-unreachable-code-fallthrough']

ext_modules = []

extensions = cythonize(
    ext_modules,
    include_path=['src', 'src/hilbert_maass'] + LIBRARY_DIRS + [SAGE_LIB],
    compiler_directives={
        'embedsignature': True,
        'language_level': '3',
    },
    gdb_debug=gdb_debug,
)

setuptools.setup(
    ext_modules=extensions,
    create_extension=create_extension,
    packages=['hilbert_maass',
              'hilbert_maass.modform'],
    dependency_links=['https://github.com/fredstro/hilbertmodgroup.git/#egg=package-1.0'],
    install_requires=[
        'hilbert-modular-group',
    ],
    package_data={
        "": ["*.pxd"]
    }
)
