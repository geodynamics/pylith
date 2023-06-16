# Running C++ unit tests and MMS tests

## Running C++ unit tests

The C++ unit tests (`libtests`) are implemented using `CppUnit` and a common test driver, `pylith::testing::TestDriver` in `tests/src/driver_cppunit.cc`.
`TestDriver` provides support for command line arguments to control the tests run and set PETSc options.
All of the C++ unit test or MMS test executables support the following command line arguments:

* **`--help`** Show help for command line arguments.
* **`--list`** List tests run by the executable.
* **`--tests=TESTS`** Run subset of the tests. `TESTS` is a comma separated list of tests.
* **`--petsc VALUE=ARG`** Set PETSc option `-VALUE=ARG`.
* **`--journal.info=NAME`** Activate Pythia info journal for `NAME`.
* **`--journal.debug=NAME`** Activate Pythia debug journal for `NAME`.
* **`--journal.warning=NAME`** Activate Pythia warning journal for `NAME`.

```{code-block} console
---
caption: Examples of using command line arguments in running C++ unit tests.
---
$ cd tests/libtests/problems

# List tests
$ ./test_problems --list

# Run all TestObserversPhysics tests.
$ ./test_problems --tests=pylith::problems::TestObserversPhysics

# Run TestObserversPhysics testVerifyObservers and testNotifyObservers tests.
$ ./test_problems --tests=pylith::problems::TestObserversPhysics::testVerifyObservers,pylith::problems::TestObserversPhysics::testNotifyObservers

# Turn on timedependent info journal.
$ ./test_problems --journal.info=timedependent

# Turn on timedependent debug journal.
$ ./test_problems --journal.debug=timedependent
```

# Running MMS tests

The Method of Manufactured Solutions (MMS) tests are implemented using `Catch2` and a common test driver, `pylith::testing::TestDriver` in `tests/src/driver_catch2.cc`.
`TestDriver` provides support for command line arguments to control the tests run and set PETSc options.
All of the C++ unit test or MMS test executables support the following command line arguments:

* **`--help`** Show help for command line arguments.
* **`--list-tests`** List tests run by the executable.
* **`--petsc VALUE=ARG`** Set PETSc option `-VALUE=ARG`.
* **`--journal.info=NAME`** Activate Pythia info journal for `NAME`.
* **`--journal.debug=NAME`** Activate Pythia debug journal for `NAME`.
* **`--journal.warning=NAME`** Activate Pythia warning journal for `NAME`.

```{code-block} console
---
caption: Examples of using command line arguments in running C++ and MMS tests.
---
$ cd tests/mmstests/linearelasticity/nofaults-2d

# List tests
$ ./mmstest_linearelasticity_nofaults2d --list-tests

# Run all UniformStrain2D tests.
$ ./mmstest_linearelasticity_nofaults2d [UniformStrain2D]

# Run all UniformStrain2D residual tests
$ ./mmstest_linearelasticity_nofaults2d [UniformStrain2D][testResidual]

# Turn on timedependent info journal.
$ ./mmstest_linearelasticity_nofaults2d --journal.info=timedependent

# Turn on timedependent debug journal.
$ ./mmstest_linearelasticity_nofaults2d --journal.debug=timedependent
```

## Using the debugger

The executables in the build directory are shell script wrappers created by `libtool`.
The underlying binary executables are in the `.libs` directory.
When using the debugger, pass the binary executable to the debugger.
For example, `gdb .libs/test_problems`.

## Using valgrind

```{code-block} console
# Run valgrind on the test_problems executable
$ valgrind --log-file=valgrind_problems.log \
  --suppressions=$PYLITH_DIR/share/valgrind-python.supp .libs/test_problems
```
