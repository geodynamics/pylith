(sec-user-run-pylith-troubleshooting)=
# Troubleshooting

## Tips and Hints For Running PyLith

* Examine the examples for a problem similar to the one you want to run and dissect it in detail.
* Start with a uniform-resolution coarse mesh to debug the problem setup.  Increase the resolution as necessary to resolve the solution fields of interest (resolving stresses/strains may require a higher resolution than that for resolving displacements).
* Merge materials using the same material model. This will result in only one VTK or HDF5 file for each material model rather than several files.
* The rate of convergence in quasistatic (implicit) problems can sometimes be improved by renumbering the vertices in the finite-element mesh to reduce the bandwidth of the sparse matrix. PyLith can use the reverse Cuthill-McKee algorithm to reorder the vertices and cells.
* If you encounter errors or warnings, run `pylith_dumpparameters` or use the `--help`, `--help-components`, or `--help-properties` command-line arguments when running PyLith to check the parameters to make sure PyLith is using the parameters you intended.
* Use the `--petsc.log_summary`, `--petsc.ksp_monitor`, `--petsc.ksp_view`, `--petsc.ksp_converged_reason`, and `--petsc.snes_converged_reason` command-line arguments (or set them in a parameter file) to view PyLith performance and monitor the convergence.
* Turn on the journals (see the examples) to monitor the progress of the code.

## Troubleshooting

Consult the PyLith category in the [CIG community forum](https://community.geodynamics.org) to see if someone else encountered a similar issue.

## Common Error Messages

### Import Error and Missing Library

```{code-block} bash
ImportError: liblapack.so.2: cannot open shared object file: No such file or directory
```

PyLith cannot find one of the libraries.
You need to set up your environment variables (e.g., `PATH`, `PYTHONPATH`, and `LD_LIBRARY_PATH`) to match your installation.
If you are using the PyLith binary on Linux or macOS, run the command `source setup.sh` in the directory where you unpacked the distribution.
This will set up your environment variables for you.
If you are building PyLith from source, please consult the instructions for building from source.

### Unrecognized Property 'p4wd'

```{code-block} bash
-- pyre.inventory(error) } \\
-- p4wd <- 'true' } \\
-- unrecognized property 'p4wd' } \\
>> command line:: } \\
-- pyre.inventory(error) } \\
-- p4pg <- 'true' } \\
-- unrecognized property ' p4pg'}
```

Verify that the `mpirun` command included in the PyLith package is the first one on your `PATH` by running the command `which mpirun`.
If it is not, adjust your `PATH` environment variable accordingly.

### Detected zero pivor in LU factorization

```{code-block} bash
-- Solving equations.
[0] PETSC ERROR: ----------------
Error Message -------------------------------
[0] PETSC ERROR: Detected zero pivot in LU factorization
see http://www.mcs.anl.gov/petsc/petsc-as/documentation/faq.html\#ZeroPivot!
```

This usually occurs when the null space of the system Jacobian is nonzero, such as the case of a problem without Dirichlet boundary conditions on any boundary.
If this arises when using the split fields and algebraic multigrid preconditioning, and no additional Dirichlet boundary conditions are desired, then the workaround is to revert to using the Additive Schwarz preconditioning without split fields as discussed in {ref}`sec-user-run-pylith-petsc-options`.

### Bus Error

This often indicates that PyLith is using incompatible versions of libraries.
This can result from changing your environment variables after configuring or installing PyLith (when building from source) or from errors in setting the environment variables `PATH`, `LD_LIBRARY_PATH`, and `PYTHONPATH`).
If the former case, simply reconfigure and rebuild PyLith.
In the latter case, check your environment variables (order matters!) to make sure PyLith finds the desired directories before system directories.

### Segmentation Fault

A segmentation fault usually results from an invalid read/write to memory.
It might be caused by an error that wasn't trapped or a bug in the code.
Please report these cases so that we can fix these problems (either trap the error and provide the user with an informative error message, or fix the bug).
If this occurs with any of the problems distributed with PyLith, simply submit a bug report (see {ref}`sec-getting-help`) indicating which problem you ran and your platform.

:::{important}
PETSc will often report errors as semgentation faults even if the underlying problem is not an invalid read/write.
If you see PETSc reporting a segmentation fault, examine the output carefully for other error messages that indicate the real issue is something else.
:::

If the crash occurs for a problem you created, it is a great help if you can try to reproduce the crash with a very simple problem (e.g., adjust the boundary conditions or other parameters of one of the examples to reproduce the segmentation fault).
Submit a bug report along with log files showing the backtrace from a debugger (e.g., gdb) and the valgrind log file (only available on Linux platforms).
You can generate a backtrace using the debugger by using the `--petsc.start_in_debugger` command-line argument:

```{code-block} console
$ pylith [..args..] --petsc.start_in_debugger
(gdb) continue
(gdb) backtrace
```

To use valgrind to detect the memory error, first go to your working directory and run the problem with `--launcher.dry`:

```{code-block} console
$ pylith [..args..] --launcher.dry
```

Instead of actually running the problem, this causes PyLith to dump the mpirun/mpiexec command it will execute.
Copy and paste this command into your shell so you can run it directly.
Insert the full path to valgrind before the full path to mpinemesis and tell valgrind to use a log file:

```{code-block} console
$ mpirun /path/to/valgrind --log-file=valgrind-log /path/to/mpinemesis --pyre-start
  [..lots of junk..]
```

% End of file
