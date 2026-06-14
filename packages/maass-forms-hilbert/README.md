# maass-forms-hilbert
Algorithms for Hilbert Maass forms

# Installing

The recommended install pulls in the database extras:

```bash
pip install 'maass_forms_hilbert[db]'
```

Plain `pip install maass_forms_hilbert` works for the math-only surface, but
`maass_forms_hilbert.database` will raise `ModuleNotFoundError` on import until
`[db]` is added. The `[db]` extras chain to `maass_form_core[db]`, which pulls
in `mongoengine`, `comp_manager`, and `httpx2`.

From a source checkout:

> make install 

ALT

> make sdist
> sage -pip install --no-build-isolation dist/*

## Extras

- **`[db]`** — MongoDB/MongoEngine persistence layer (via `maass_form_core[db]`). Required for `maass_forms_hilbert.database`.
- **`[dev]`** — Tooling for contributors (`pytest`, `ruff`, `tox`, `pre-commit`, …).

# Testing
1. Doctests
2. Pytest
    To enable pytests you need to first run `sage -i pytest pytest_xdist`