% -*- markdown -*-
% -*- coding: utf-8 -*-
%
% This code is part of PyLith, developed through the Computational Infrastructure
% for Geodynamics (https://github.com/geodynamics/pylith).
%
% Copyright (c) 2010-2026 The PyLith authors and contributors
% SPDX-License-Identifier: MIT

# Coding Style

There are a number of standard coding styles for programming languages, notably PEP8 for Python. For PyLith, we try to be consistent in naming conventions across Python and C++ while following a subset of the conventions used in PETSc and PEP8.

:::{important}
We use 4 spaces for indentation. Configure your editor to use spaces instead of tabs.
:::

:::{tip}
Run the Python script `etc/ci/conventions/check.py` to verify your files match basic style conventions BEFORE submitting a pull request.
The continuous integration pipeline will also run this script and fail if code does not match the basic style conventions.
:::



## General guidelines

- Naming conventions
  - Use self-documenting names.
  - Avoid single letter variables. Choose meaningful names that will be found via searches across single or multiple files (e.g., grep).
  - Class names are generally nouns and methods are verbs.
  - Class names use upper camel case, e.g., `TimeDependent`.
  - Public method names use camel case, e.g., `computeRHSResidual()`.
  - Protected and private method names use camel case preceded by an underscore, e.g., `_setFEKernelsRHSResidual()`.
  - In C++ data members are private and use camel case preceded by an underscore, e.g., `_gravityField`.
  - In Python data members are public and use camel case, e.g., `self.gravityField`.
  - Local variables use camel case, e.g., `numIntegrators`.
- Comments
  - Include the standard banner at the very beginning of every file.
  - For every class method and function, describe its function using a narrative that mentions each argument (marked by braces).
  - Document nontrivial algorithms and any assumptions.
- Error checking
  - PyLith should never crash without an error message.
  - All user errors should be trapped as early as possible and reported with an informative error message by throwing an appropriate exception. If possible, suggest ways to correct the error.
  - Messages for internal errors should indicate the location in the code where the error was trapped.
  - All pointers should be checked for `nullptr` values before use. Usually we use `assert()` to do this.
  - Check the return values for all calls to functions in external libraries. For PETSc functions, we use `PylithCallPetsc(PetscFunction(args));`.
- Commit messages
  - Start with an identifier, e.g., `lib`, `ci`, `libtest`, `mms`, `dmplex`, `petsc`, following by a colon and a descriptive message describing the change.
  - Multiple identifiers are separated by commas.
- Testing
  - All C++ methods should be covered by unit tests.
  - All governing equations should be covered by Method of Manufactured solution tests.
  - All functionality should be covered by full-scale tests.

## Formatting source files

We use `black --line-length=120` and `clang-format` (`.clang-format` in the top-level directory is read automatically) to format Python and C/C++ source files, respectively.

:::{tip}
We highly recommend using an integrated development environment, such as Visual Studio Code, that allow `black` and `clang-format` to automatically format all C/C++ and Python source files.
:::

## Error Checking

Our philosophy is that PyLith should never crash without an error message.
If it encounters a fatal error, then it should generate an appropriate error message and abort.
In C++ we throw specialized exceptions (`libsrc/pylith/exceptions/Exceptions.hh`) based on the type of error.

Additional protections against crashing include: using asserts to verify pointers are non-null before using them and using the `PylithCallPetsc()` macro to check the return value for *every* call to a PETSc function.

```{code-block} c++
---
caption: Example of using `assert()`
emphasize-lines: 1
---

assert(_solution); // Verify _solution is not NULL.

// Initialize integrators.
const size_t numIntegrators = _integrators.size();
for (size_t i = 0; i < numIntegrators; ++i) {
    assert(_integrators[i]); // Verify _integrators[i] is not NULL.
    _integrators[i]->initialize(*_solution);
} // for  
```

```{code-block} c++
---
caption: Example of using `PylithCallPetsc()` macro.
---
PylithCallPetsc(TSGetTimeStep(ts, &dt));
```

:::{tip}
When we build the code for production runs, we usually configure with `CPPFLAGS=-DNDEBUG` to remove `assert()` calls.
:::

In combination with the above procedures, we also make use of the Pyre journals to display warnings and errors to facilitate debugging.
The journals provide the file name and line number along with the message.
By default, Pyre journals for errors and warnings are turned on and those for debugging are turned off.
The header file `utils/journals.hh` defines several useful macros for reporting warnings and errors.

```{code-block} c++
---
caption: Example of using Pyre journals and standard exceptions.
emphasize-lines: 9, 12-13
---
switch (bitUse) {
case 0x1:
    _bcKernel = pylith::fekernels::TimeDependentFn::initial_scalar;
    break;
case 0x2:
    _bcKernel = pylith::fekernels::TimeDependentFn::rate_scalar;
    break;
case 0x0:
    PYLITH_COMPONENT_WARNING("Dirichlet BC provides no constraints.");
    break;
default:
    PYLITH_ERROR(pylith::exceptions::InternalLogicError,
        "Unknown combination of flags for Dirichlet BC terms "
        << "(useInitial="<<_useInitial<<", useRate="<<_useRate<<").");
} // switch
```

## C/C++ style

### Object Declaration Files

C++ object declaration (header) files use the `.hh` suffix.
C header files use the `.h` suffix.

### Object Implementation Files

C++ object implementation files use the `.cc` suffix.
Inline implementation files use the `.icc` suffix and are included from the definition (header) files.
C implementation files use the `.c` suffix.

To facilitate debugging and error messages, we use the following
macros:

`PYLITH_METHOD_BEGIN`
: This macro allows line numbers of source files to be included in PETSc error messages. Use this macro at the beginning of all methods using any PETSc routines as well as most other methods. We don't use this macro in destructors because many of them are called *after* `PetscFinalize`. We also do not use this macro in trivial or inline methods that do not call any PETSc routines.

`PYLITH_METHOD_END`
: Use the macro at the end of all methods that begin with `PYLITH_METHOD_BEGIN` and return void.

`PYLITH_RETURN_END`
: Use this macro at the end of all methods that begin with `PYLITH_METHOD_BEGIN` and return non-void values.

`PylithCallPetsc()`
: Use this macro in *every* call to a PETSc function to check the return value.

`PYLITH_JOURNAL_DEBUG`
: Use this macro immediately after `PYLITH_METHOD_BEGIN` in methods of all objects that inherit from `GenericComponent`.

`PYLITH_COMPONENT_DEBUG`
: Use this macro immediately after `PYLITH_METHOD_BEGIN` in methods of all objects that inherit from `PyreComponent`.
Non-abstract classes should call `PyreComponent::setName()` in the constructor.
We recommend using a static data member for the name with the lowercase name matching the Pyre component, e.g., "timedependent" for the C++ `TimeDependent` object.

## Python style

Python files use the `.py` suffix.
All classes, methods, and functions should be documented using docstrings.

% End of file
