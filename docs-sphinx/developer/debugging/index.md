# Debugging

## C++ unit tests and MMS tests

The C++ unit tests (`libtests`) and Method of Manufactured Solutions (MMS) tests are implemented using `CppUnit` and a common test driver, `pylith::testing::TestDriver`. `TestDriver` provides support for command line arguments to control the tests run and set PETSc options. All of the C++ unit test or MMS test executables support the following command line arguments:

* **`--help`** Show help for command line arguments.
* **`--list`** List tests run by the executable.
* **`--tests=TESTS`** Run subset of the tests. `TESTS` is a comma separated list of tests.
* **`--petsc VALUE=ARG`** Set PETSc option `-VALUE=ARG`.

Example
```

```

:::{admonition} Using the debugger with C++ unit tests and MMS tests
:class: important

The executables in the build directory are shell script wrappers created by `libtool`. The underlying binary executables are in the `.libs` directory. When using the debugger, pass the binary executable to the debugger. For example, `gdb .libs/test_feassemble`.
:::
