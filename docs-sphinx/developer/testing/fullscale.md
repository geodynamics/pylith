# Full-scale tests

:::{error}
TODO Add overview of full-scale testing.
:::

## Using the debugger

To start the `gdb` debugger when running the PyLith application, simply add the command line argument `--petsc.start_in_debugger`.
To use an alternative debugger, such as `lldb`, append the name of the debugger executable, for example `--petsc.start_in_debugger=lldb`.
By default, PETSc will try to start the debugger in an xterm.
To use an alternative terminal program, use the command line argument `--petsc.debug_terminal=TERMINAL`.
For example for the GNOME terminal, use `--petsc.debug_terminal="gnome-terminal -x"`.

## Using valgrind

Running valgrind on a the PyLith application, including full-scale tests, requires using an additional argument (`--trace-children=yes`), because the PyLith application forks a subprocess to do the computation.

```{code-block} console
# Run valgrind on the PyLith executable
$ valgrind --trace-children=yes \
  --suppressions=$PYLITH_DIR/share/valgrind-python.supp pylith step01.cfg
```
