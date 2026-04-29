# maass-forms

A mathematical software ecosystem for Maass forms computation, consisting of domain-specific packages and shared infrastructure built on SageMath.

## Packages

| Package | Description |
|---------|-------------|
| `maass-form-core` | Shared infrastructure (database, coefficients, utilities) |
| `maass-forms-hilbert` | Hilbert Maass forms over number fields |
| `maass-forms-klein` | Kleinian Maass forms on hyperbolic 3-manifolds |

## Installation

```bash
# Install all packages (core first, then domain packages)
make install

# Or install individually
cd packages/maass-forms-hilbert && make install
```

## Versioning

Each package is versioned independently using [setuptools-scm](https://github.com/pypa/setuptools-scm) with **prefixed git tags**:

| Package | Tag prefix | Example tag |
|---------|-----------|-------------|
| `maass-form-core` | `core/v` | `core/v0.1.0` |
| `maass-forms-hilbert` | `hilbert/v` | `hilbert/v0.2.0` |
| `maass-forms-klein` | `klein/v` | `klein/v0.1.3` |

### Creating a release

Tag the commit with the appropriate prefix:

```bash
# Release core v0.1.0
git tag -a core/v0.1.0 -m "maass-form-core 0.1.0"

# Release hilbert v0.2.0
git tag -a hilbert/v0.2.0 -m "maass-forms-hilbert 0.2.0"

# Release klein v0.1.3
git tag -a klein/v0.1.3 -m "maass-forms-klein 0.1.3"

git push --tags
```

Multiple packages can be released from the same commit if needed.

### How it works

Each package's `pyproject.toml` configures `setuptools_scm` with:
- `git_describe_command` using `--match <prefix>/v*` to find only its own tags
- `tag_regex` to strip the prefix and extract the version number
- `fallback_version = "0.0.0"` for builds before the first tag

Between tags, setuptools-scm generates development versions automatically (e.g. `0.2.0.dev3+gabcdef`).

## Development

```bash
make test      # Run tests across all packages
make lint      # Lint all packages with ruff
make format    # Format all packages with ruff
make clean     # Clean build artifacts
```

See the per-package `CLAUDE.md` files for detailed development instructions.

## License

GPL-3.0-or-later
