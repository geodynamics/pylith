# PyLith

PyLith is an application for parallel numerical modeling of crustal deformation with a focus on earthquake faulting.
The top-level code is written in Python with underlying code in C++ with Python bindings written using SWIG.
This file must be self-contained.
Do not rely on linked Markdown files being read automatically.
The essential repo guidance is embedded below.

## Project Layout

- `ci-config/` -  Configuration files and scripts used in continuous integration
- `developer/` - Configuration files and scripts used by developers
- `docker/` - Files for building Docker images
- `examples/` - Input files for running PyLith examples
- `libsrc/` - C++ source files
- `modulesrc/` - Python bindings
- `pylith` - Python source code
- `docs/` - user and developer documentation
- `release-notes/` - Directory with release notes
- `share/` - Directory with files distributed with installation (generic parameter files)
- `templates/` - Directory with examples demonstrating how to extend the code
- `tests/libtests` - Unit tests for the C++ library
- `tests/mmstest` - MMS tests for the C++ library
- `tests/pytests` - Basic tests for the Python objects and SWIG interface
- `tests/fullscale` - Full-scale tests that run PyLith and verify the output
- `benchmarks/` - Solver and performance benchmarks
- `playpen/` - Directory containing old, experimental tests and code

## Core Working Rules

- Preserve PyLith style and naming conventions.
- Keep edits minimal and local to the requested change.
- Match existing patterns in the directory you are modifying before introducing a new one.
- Avoid unnecessary code duplication. Prefer reusing or extending nearby logic when it keeps behavior clear and local.
- Do not add speculative abstractions or broad refactors unless explicitly requested.
- If you modify source, check whether test blocks, expected output files, or documentation need corresponding updates.

## PyLith Naming And API Conventions

- Enum constants and macros are uppercase with underscores, for example `PYLITH_METHOD_BEGIN`.
- Class data members begin with underscores such as `_problem`.
- Class private and protected members begin with underscores such as `getDM()`.

## PETSc Data Type Rules

- Use `PylithInt` for most indices and array lengths.
- Use `size_t` for sizes.
- Do not silence narrowing warnings with blind casts.
- Prefer PETSc MPI wrappers that accept PETSc count types when large counts may be involved.

## C++ Coding Style

- Formatting is controlled by `.clang-format`. Use `make clangformat` when needed.
- Group local variables by type with one variable per line.
- Initialize local variables in the declaration when practical.
- Begin every class method with `PYLITH_METHOD_BEGIN;` and end with `PYLITH_METHOD_END;` for void functions and `PYLITH_METHOD_RETURN(value);` for functions returning values.
- For `PetscErrorCode` functions, return `PetscFunctionReturn(PETSC_SUCCESS)` on success.
- Wrap PETSc calls with `PyLithCallPetsc(...)`.
- Use the macros in `journals.hh` for printing to stdout for informational messages, debugging messages, and warnings. Use the `PYLITH_COMPONENT_*` version for `PyreComponent` objects. Use the simpler version such as `PYLITH_WARNING` for other obejcts. 
- Do not leave commented-out code or dead `#ifdef` blocks in source files.
- Use `/* ... */` for multiline comments and `// ...` for short single-line comments.
- Decorate multiline comments with leading `*` on each line.
- Always append `()` to function names when mentioning them in comments, for example `getValue()`.
- Use correct grammar and spelling in comments and messages.

## Error Handling And PETSc Idioms

- Most PETSc functions return `PetscErrorCode`.
- Use `assert()` to check pointers before using them.
- Reuse existing PyLith utility routines and macros before adding custom helpers.

## Testing Requirements

- If behavior changes, update the corresponding tests.
- Data files for tests live in `data` under the test directory.
- Use Catch2 as the C++ test framework with tests implemented in a class for each C++ class being tested.
- Keep tests targeted. Add or update the narrowest test that proves the behavior you changed.

## Build And Test Commands

- `make` - Build libraries and SWIG modules
- `make install` - Install libraries, SWIG modules, and Python code
- `make check` - run tests

## Merge Request Expectations

- All changes are expected to arrive through GitHub pull requests.
- Keep diffs reviewable and focused.
- Before concluding work, consider whether formatting, source-style checks, and at least one relevant test should be run.
- If you cannot run the appropriate verification in the current environment, say so explicitly.

## Practical Agent Guidance

- Read nearby code before editing so new code matches local conventions.
- When touching PyLith C++ code, check for consistent use of `PyLithCallPetsc`, `PYLITH_METHOD_BEGIN`, naming, and test coverage.
- When touching tutorials or tests, inspect neighboring files for the expected structure and output-file conventions.
- When touching public interfaces, check whether headers, docs, and examples need updates.
- Prefer citing exact file paths and commands in your responses.

## Anti-Patterns (MUST avoid when writing or reviewing)

## Key References

External links are for human convenience only; do not assume linked Markdown files will be ingested automatically.

- Release docs: https://pylith.readthedocs.io
- GitLab project: https://github.com/geodynamics/pylith
