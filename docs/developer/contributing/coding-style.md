# Coding Style

There are a number of standard coding styles for programming languages, notably PEP8 for Python. For PyLith, we try to be consistent in naming conventions across Python and C++ while following a subset of the conventions used in PETSc and PEP8 with documentation styles consistent with Doxygen.

:::{important}
We use 4 spaces for indentation. Configure your editor to use spaces instead of tabs.
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
  - Include standard banner at the very beginning of every file.
  - For every class method, describe its function and include a description for each argument. For Python this is done in the docstring of the method, and for C++ this is done in a doxygen style comment immediately before the method declaration in the header file.
  - Document nontrivial algorithms and any assumptions.
- Error checking
  - PyLith should never crash without an error message.
  - All user errors should be trapped as early as possible and reported with an informative error message by throwing an appropriate exception. If possible, suggest ways to correct the error.
  - Messages for internal errors should indicate the location in the code where the error was trapped.
  - All pointers should be checked for `NULL` values before use. Usually we use `assert()` to do this.
  - Check the return values for all calls to functions in external libraries. For PETSc functions, we use `PYLITH_CHECK_ERROR(returnValue)`.
- Testing
  - All C++ methods should be covered by unit tests.
  - All governing equations should be covered by Method of Manufactured solution tests.
  - All functionality should be covered by full-scale tests.

## Formatting source files

We use `autopep8` and `uncrustify` to format Python and C/C++ source files, respectively.
The corresponding configuration files are `autopep8.cfg` and `uncrustify.cfg` in the `developer` directory.
The Python script `developer/format_source.py` is a handy utility for calling `autopep8` and `uncrustify` with the appropriate arguments and formatting multiple files.

:::{tip}
We highly recommend using an integrated development environment, such as Visual Studio Code, that allow `uncrustify` and `autopep8` to automatically format all C/C++ and Python source files.
:::

```{code-block} bash
---
caption: Formatting Python and C++ source code using `format\_source.py`. `autopep8` and `uncrustify` must be in the current path.
---
# Format a C++ file.
developer/format_source.py --cplusplus=libsrc/pylith/materials/Material.cc

# Format all C++ files in the 'libsrc/pylith/materials' directory.
developer/format_source.py --cplusplus=libsrc/pylith/materials/*.cc

# Format a Python file.
developer/format_source.py --python=pylith/materials/Material.py

# Format all Python files in the 'pylith/materials' directory.
developer/format_source.py --python=pylith/materials/*.py
```

## Error Checking

Our philosophy is that PyLith should never crash without an error message.
If it encounters a fatal error, then it should generate an appropriate error message and abort.
In C++ we throw `std::runtime_error` exceptions for errors resulting from user input and `std::logic_error` exceptions for internal inconsistencies or logic errors.
In Python we use standard exception objects.

Additional protections against crashing include: using asserts to verify pointers are non-null before using them and using the `PYLITH_CHECK_ERROR` macro to check the return value after *every* call to a PETSc function.

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

:::{tip}
When we build the code for production runs, we usually configure with `CPPFLAGS=-DNDEBUG` to remove `assert()` calls.
:::

```{code-block} c++
---
caption: Example of using `PYLITH_CHECK_ERROR` macro.
---
PetscErrorCode err = TSGetTimeStep(ts, &dt);PYLITH_CHECK_ERROR(err);
```

In combination with the above procedures, we also make use of the Pyre journals to display warnings and errors to facilitate debugging.
The journals provide the file name and line number along with the message.
By default, Pyre journals for errors are turned on and those for warnings and debugging are turned off.
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
    PYLITH_COMPONENT_LOGICERROR("Unknown combination of flags for Dirichlet BC terms "
        << "(useInitial="<<_useInitial<<", useRate="<<_useRate<<").");
} // switch
```

## C/C++ style

### Object Declaration Files

C++ object declaration (header) files use the `.hh` suffix.
C header files use the `.h` suffix.

:::{important}
*All* declarations of class methods should include a description of what the method does and a description of each argument and the return value if it is not void.
:::

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

`PYLITH_CHECK_ERROR`
: Use this macro after *every* call to a PETSc function to check the return value.

`PYLITH_JOURNAL_DEBUG`
: Use this macro immediately after `PYLITH_METHOD_BEGIN` in methods of all objects that inherit from `GenericComponent`.

`PYLITH_COMPONENT_DEBUG`
: Use this macro immediately after `PYLITH_METHOD_BEGIN` in methods of all objects that inherit from `PyreComponent`.
Non-abstract classes should call `PyreComponent::setName()` in the constructor.
We recommend using a static data member for the name with the lowercase name matching the Pyre component, e.g., "timedependent" for the C++ `TimeDependent` object.

## Python style

Python files use the `.py` suffix.
All classes, methods, and functions should be documented using docstrings.
