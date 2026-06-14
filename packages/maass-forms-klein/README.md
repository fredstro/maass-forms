# maass-forms-klein
Algorithms for Maass forms on knot complements


## Requirements
- SageMath v10.3+ (https://www.sagemath.org/)
- Snappy 3.0.3  (run `sage -i snappy` to install it)


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

The following commands are run in an isolated docker container 
and requires docker to be installed and running:
1. `docker` -- build a docker container with the tag `hilbertmodgroup-{GIT_BRANCH}`
2. `docker-rebuild` -- rebuild the docker container without cache
3. `docker-test` -- run SageMath's doctests in the docker container
4. `docker-examples` -- run a Jupyter notebook with the SageMath kernel initialised at the `/examples` directory 
  and exposing the notebook at http://127.0.0.1:8888. The port used can be modified by 
5. `docker-tox` -- run tox with all environments: `doctest`, `coverage`, `pycodestyle`, `relint`, `codespell`. 
6. `docker-shell` -- run a shell in a docker container
7. `docker-sage` -- run a sage interactive shell in a docker container

The following command-line parameters are available 
- `NBPORT` -- set the port of the notebook for `examples` and `docker-examples`  (default is 8888)
- `TOX_ARGS` -- can be used to select one or more of the tox environments (default is all)
- `REMOTE_SRC` -- set to 0 if you want to use the local source instead of pulling from gitHub (default 1)
- `GIT_BRANCH` -- the branch to pull from gitHub (used if REMOTE_SRC=1)

### Example usage
Run tox coverage on the branch `main` from gitHub:

`make docker-tox REMOTE_SRC=1 GIT_BRANCH=main TOX_ARGS=coverage`

Run doctests on the local source with local version of sage:

`make tox TOX_ARGS=doctest`

Run relint on the local source with docker version of sage:

`make docker-tox REMOTE_SRC=0 TOX_ARGS=relint`

Run a jupyter notebook server with a sage kernel that has `maass_forms_klein` installed
and uses the local source files (note that the )

`make docker-examples REMOTE_SRC=0`

# Installation

The recommended install pulls in the database extras:

```bash
pip install 'maass_forms_klein[db]'
```

Plain `pip install maass_forms_klein` works for the math-only modules, but
`maass_forms_klein.database` will raise `ModuleNotFoundError` on import until
`[db]` is added. The `[db]` extras chain to `maass_form_core[db]`, which pulls
in `mongoengine`, `comp_manager`, and `httpx2`.

From a source checkout, run `make install`.

## Extras

- **`[db]`** — MongoDB/MongoEngine persistence layer (via `maass_form_core[db]`). Required for `maass_forms_klein.database`.
- **`[dev]`** — Tooling for contributors (`pytest`, `ruff`, `tox`, `pre-commit`, …).
- **`[notebook]`** — Jupyter notebook server for examples.

# Testing
Run `make test`

