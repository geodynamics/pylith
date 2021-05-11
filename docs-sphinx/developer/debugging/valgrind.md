# Using valgrind

Valgrind is a useful tool for finding memory leaks, use of uninitialized variables, and invalid reads and writes to memory.
When running valgrind there are three very useful command line arguments:

* **`--log-filename=FILENAME`** Send output to FILENAME. This does not work when running the PyLith
application because each new process wipes out the log file.
* **`----suppressions=FILE`** Omit errors matching given patterns when reporting errors. Valgrind often reports lots of errors arising from the way OpenMPI and Python handle memory allocation and deallocation. We usually use the Python suppression file `share/valgrind-python.supp` when running valgrind.
* **`--trace-children=yes`** Continue tracing errors in subprocesses. This is important when running valgrind on the PyLith executable, as the actual computation is done in a forked process.

## C++ and MMS tests

```{code-block} console
# Run valgrind on the test_problems executable
$ valgrind --log-file=valgrind_problems.log \
  --suppressions=$PYLITH_DIR/share/valgrind-python.supp .libs/test_problems
```

## Full-scale tests

Running valgrind on a the PyLith application, including full-scale tests, requires using an additional argument (`--trace-children=yes`), because the PyLith application forks a subprocess to do the computation.

```{code-block} console
# Run valgrind on the PyLith executable
$ valgrind --trace-children=yes \
  --suppressions=$PYLITH_DIR/share/valgrind-python.supp pylith step01.cfg
```
