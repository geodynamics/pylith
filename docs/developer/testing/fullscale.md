# Full-scale tests

The full-scale tests are constructed in Python, making use of the standard Python `unittest` module.
A full-scale test involves running a PyLith simulation and using Python code to check output against a known solution.
Whenever possible, we simulate simple boundary value problems with analytical solutions.
We run some full-scale tests in parallel on a small number of processes.

The first step is to construct a PyLith simulation (mesh, parameter files, and spatial databases) that solves the desired boundary value problem with a known solution.
In constructing the simulation, we create Python scripts that generate the spatial database files.
This makes it easy to synchronize any changes between the problem specification and the known solution.

Once we have a PyLith simulation that appears to run correctly, we construct a Python object that inherits from `pylith.testing.FullTestApp`.
See `tests/fullscale/linearelasticity/TestAxialDisp.py` for an example.
This Python object sets up the a PyLith simulation, runs the simulation, and defines methods that check fields in the output files.
In most cases we solve the boundary value problems in 2D with both quadrilateral and triangular meshes and in 3D with both tetrahedral and hexahedral meshes.

When multiple full-scale tests use the same meshes, we place them in the same directory and use dictionaries in a `meshes.py` file to define the mesh information.

## Command line arguments

The `pylith.testing.FullTestApp` Python application has two optional command line arguments:

* **`--verbose` (int)** Set verbosity level. Level > 0 will display which fields are being checked.
* **`--skip-pylith-run` (bool)** Skip running the PyLith simulation and only check the existing output files against the expected results. This argument is useful when debugging the checks.


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

## Example

:::{admonition} TODO
:class: error

Walk through an example of a full-scale test.
:::
