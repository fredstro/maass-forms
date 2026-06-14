# maass-form-core

Shared infrastructure library for the Maass forms computation ecosystem.

## Purpose

This package extracts common patterns from `maass-forms-hilbert` and `maass-forms-klein` to eliminate ~70% code duplication while maintaining mathematical domain separation.

## Installation

The recommended install pulls in the database extras:

```bash
pip install 'maass_form_core[db]'
```

Plain `pip install maass_form_core` works for the math-only surface, but `maass_form_core.database` will raise `ModuleNotFoundError` on import until `[db]` is added. The domain packages (`maass-forms-hilbert`, `maass-forms-klein`) use the database layer, so install with `[db]` unless you have a specific reason not to.

### Extras

- **`[db]`** — MongoDB/MongoEngine persistence layer (`mongoengine`, `comp_manager`, `httpx2`). Required for `maass_form_core.database`.
- **`[dev]`** — Tooling for contributors (`pytest`, `ruff`, `tox`, `pre-commit`, …).

## Modules

- **database/** - Common MongoDB/MongoEngine patterns, base document classes, QuerySets *(requires `[db]`)*
- **coefficients/** - Coefficient indexing, storage, and serialization
- **spaces/** - Abstract mathematical space and form element base classes
- **functions/** - Shared special functions (Bessel functions, etc.)
- **utils/** - Type definitions (`Real_t`, `Complex_t`), validation, JSON converters
- **testing/** - Shared test fixtures and assertion helpers

## Status

This package is under active development. See `IMPLEMENTATION_PLAN_MAASS_FORMS_CORE.md` in the repository root for the migration roadmap.
