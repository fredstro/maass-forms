# knot-maass
Algorithms for Maass forms on knot complements


## Requirements
- SageMath v9.6+ (https://www.sagemath.org/)

## Development and testing

The make file `Makefile` contains a number of useful commands that you can run using 
```console
$ make <command>
```
The following commands are run in your local SagMath environment:
1. `build` -- builds the package in place (sometimes useful for development).
2. `sdist` -- create a source distribution in /sdist (can be installed using `sage -pip install sdist/<dist name>`)
3. `install` -- build and install the package in the currently active sage environment
4. `clean` -- remove all build and temporary files
5. `test` -- run sage's doctests (same as `sage -t src/*`)
6. `examples` -- run a Jupyter notebook with the SageMath kernel initialised at the `/examples` directory.
7. `tox` -- run `sage -tox` with all environments: `doctest`, `coverage`, `pycodestyle`, `relint`, `codespell`
   Note: If your local SageMath installation does not contain tox this will run `sage -pip install tox`.

# Installation
Run `make install`

# Testing
Run `make test`

