# maass-form-core

Shared infrastructure library for the Maass forms computation ecosystem.

## Purpose

This package extracts common patterns from `maass-forms-hilbert` and `maass-forms-klein` to eliminate ~70% code duplication while maintaining mathematical domain separation.

## Modules

- **database/** - Common MongoDB/MongoEngine patterns, base document classes, QuerySets
- **coefficients/** - Coefficient indexing, storage, and serialization
- **spaces/** - Abstract mathematical space and form element base classes
- **functions/** - Shared special functions (Bessel functions, etc.)
- **utils/** - Type definitions (`Real_t`, `Complex_t`), validation, JSON converters
- **testing/** - Shared test fixtures and assertion helpers

## Status

This package is under active development. See `IMPLEMENTATION_PLAN_MAASS_FORMS_CORE.md` in the repository root for the migration roadmap.
