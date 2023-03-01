import os
import setuptools
from setuptools.extension import Extension
from Cython.Build import cythonize
# Check if we are currently in a SageMath environment.
SAGE_LOCAL = os.getenv('SAGE_LOCAL')
if not SAGE_LOCAL:
    raise ValueError("This package can only be installed inside SageMath (http://www.sagemath.org)")
# Find correct value for SAGE_LIB which is needed to compile the Cython extensions.
SAGE_LIB = os.getenv('SAGE_LIB')
if not SAGE_LIB:
    try:
        from sage.env import SAGE_LIB
    except ModuleNotFoundError:
        raise ModuleNotFoundError("To install this package you need to either specify the "
                                  "environment variable 'SAGE_LIB' or call pip with "
                                  "'--no-build-isolation'")
if not os.path.isdir(SAGE_LIB):
    raise ValueError(f"The library path {SAGE_LIB} is not a directory.")

# lib_headers = {"gmp": [os.path.join(SAGE_INC, 'gmp.h')],  # cf. #8664, #9896
import pprint
# pprint.pprint(os.environ)
print("SAGE_LIB-",SAGE_LIB)
# Extension modules using Cython
extra_compile_args = ['-Wno-unused-function',
                      '-Wno-implicit-function-declaration',
                      '-Wno-unused-variable',
                      '-Wno-deprecated-declarations',
                      '-Wno-deprecated-register']
ext_modules = [
    Extension(
        'knot_maass.functions.besselk',
        sources=[os.path.join('src/knot_maass/functions/besselk.pyx')],
        extra_compile_args=extra_compile_args,
        libraries=['gmp']
    )]

print("ext modules=",ext_modules)
setuptools.setup(
    ext_modules=cythonize(
        ext_modules,
        include_path=['src', SAGE_LIB],
        compiler_directives={
            'embedsignature': True,
            'language_level': '3',
        },
    ),
)
