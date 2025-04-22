(developer-debugging-tools)=
# Debugging tools

## VS Code integration

You can use the VS Code debugger with C++ unit tests (libtests), MMS Tests, and PyLith simulations.

### C++ unit test

1. Click on `Run and Debug` in the left sidebar of the VS Code window.
2. At the top of the panel, select the test suite to run.
3. Click on the green triangle next to the name of the selected test suite.

### MMS test

1. Click on `Testing` in the left sidebar of the VS Code window.
2. Click on the `Debug Test` icon that pops up when you mouse over the name of the test you want to run.

:::{note}
Once we migrate the C++ unit tests to Catch2, they will also show up in list of tests under `Testing`.
:::

### PyLith simulation

We provide a configuration in `launch.json` to attach the VS Code debugger to a running process, such as a PyLith simulation.

1. To attach the debugger, start the PyLith simulation while adding the command line arguments `--petsc.stop_for_debugger --petsc.pause_debugger=10` (10 corresponds to the number of seconds the program will wait for you to attach the debugger before continuing).
2. Note the process id that is printed to stdout.
3. Select `Run and Debug` in the left sidebar of the VS Code window.
4. Select `Attach to process` at the top of the panel.
5. Select the `mpinemesis` process with the process id that matches the one printed to stdout.

:::{tip}
These PETSc command line arguments used to attach the debugger are provided in `$PREFIX/share/settings/attach_debugger.cfg` where `$PREFIX` is the directory where PyLith is installed.
:::

## Debugger quick reference

Please see the `gdb` and `lldb` documentation for detailed instructions.
Here we illustrate some common basic commands.

```{code-block} console
---
caption: Debugging with gdb
---
# Set breakpoint at line 150 of Material.cc
(gdb) b Material.cc:150

# Set breakpoint at exception throw
(gdb) catch throw

# Show arguments for the current frame
(gdb) info args

# Show local variables for the current frame
(gdb) info locals

# Show the contents of a local variable: p VARIABLE
(gdb) p numFields

# Show the contents of local array: p POINTER[0]@SIZE
# Print array of 4 values pointed to by variable values
(gdb) p values[0]@4

# Print stack trace
(gdb) backtrace
```

```{code-block} console
---
caption: Debugging with lldb
---
(lldb) b Material.cc:150

# Set breakpoint at exception throw
(lldb) break set -E C++

# Show local variables
(lldb) frame variable

# Show the contents of a local variable: frame variable VARIABLE
(lldb) frame variable numFields
# Alternatively
(lldb) p numFields

# Show the contents of an array of values: parray SIZE POINTER
# Show the contents of array of 10 values pointed to by the variable values.
(lldb) parray 10 values
```

```{code-block} console
---
caption: Printing PETsc and PyLith objects (same for gdb and lldb)
---
# Print PETSc section "section"
(gdb) call PetscSectionView(section, 0)

# Print PETSc DS "ds"
(gdb) call PetscDSView(ds, 0)

# Print PyLith mesh "mesh"
# Formats:
#   "::ascii_info_detail" - ASCII to stdout
#   ":mesh.txt:ascii_info_detail" - ASCII to mesh.txt
#   ":mesh.tex:ascii_latex" - LaTeX to mesh.tex
#   "vtk:mesh.vtu:vtk_vtu" - VTK to mesh.vtu
(gdb) call mesh.view("::ascii_info_detail")

# Print PyLith field "field"
# Options values:
#   0: metadata only
#   1: Metadata and section
#   2: Metadata and PETSc vector
#   3: Metadata, section, and PETSc vector
call field.view("MY LABEL", pylith::topology::Field::ViewOptions(1))
```

## Valgrind quick reference

Valgrind is a useful tool for finding memory leaks, use of uninitialized variables, and invalid reads and writes to memory.
When running valgrind there are three very useful command line arguments:

* **`--log-filename=FILENAME`** Send output to FILENAME. This does not work when running the PyLith
application because each new process wipes out the log file.
* **`----suppressions=FILE`** Omit errors matching given patterns when reporting errors. Valgrind often reports lots of errors arising from the way OpenMPI and Python handle memory allocation and deallocation. We usually use the Python suppression file `share/valgrind-python.supp` when running valgrind.
* **`--trace-children=yes`** Continue tracing errors in subprocesses. This is important when running valgrind on the PyLith executable, as the actual computation is done in a forked process.
