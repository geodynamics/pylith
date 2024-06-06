# Step 6: Error 10

## Error Message

```{code-block} console
---
caption: Error message 10 when running Step 6.
linenos: True
emphasize-lines: 5
---
$ pylith step06_twofaults.cfg

[0]PETSC ERROR: --------------------- Error Message --------------------------------------------------------------
[0]PETSC ERROR: Error in external library
[0]PETSC ERROR: Could not find values for initiation_time at (  -21860.3  -27621) in spatial database 'Fault rupture for main fault'.
[0]PETSC ERROR: WARNING! There are unused option(s) set! Could be the program crashed before usage or a spelling mistake, etc!
[0]PETSC ERROR:   Option left: name:-ksp_atol value: 1.0e-12 source: code
[0]PETSC ERROR:   Option left: name:-ksp_converged_reason (no value) source: code
[0]PETSC ERROR:   Option left: name:-ksp_error_if_not_converged (no value) source: code
[0]PETSC ERROR:   Option left: name:-ksp_guess_pod_size value: 8 source: code
[0]PETSC ERROR:   Option left: name:-ksp_guess_type value: pod source: code
[0]PETSC ERROR:   Option left: name:-ksp_rtol value: 1.0e-12 source: code
[0]PETSC ERROR:   Option left: name:-mg_fine_pc_type value: vpbjacobi source: code
[0]PETSC ERROR:   Option left: name:-pc_type value: gamg source: code
[0]PETSC ERROR:   Option left: name:-snes_atol value: 1.0e-9 source: code
[0]PETSC ERROR:   Option left: name:-snes_converged_reason (no value) source: code
[0]PETSC ERROR:   Option left: name:-snes_error_if_not_converged (no value) source: code
[0]PETSC ERROR:   Option left: name:-snes_max_it value: 1 source: command line
[0]PETSC ERROR:   Option left: name:-snes_monitor (no value) source: code
[0]PETSC ERROR:   Option left: name:-snes_rtol value: 1.0e-12 source: code
[0]PETSC ERROR:   Option left: name:-ts_error_if_step_fails (no value) source: code
[0]PETSC ERROR:   Option left: name:-ts_monitor (no value) source: code
[0]PETSC ERROR:   Option left: name:-ts_type value: beuler source: code
[0]PETSC ERROR: See https://petsc.org/release/faq/ for trouble shooting.
[0]PETSC ERROR: Petsc Development GIT revision: v3.21.2-167-g4fed2113cae  GIT Date: 2024-05-31 10:11:14 -0400
[0]PETSC ERROR: /software/unix/py3.12-venv/pylith-debug/bin/mpinemesis on a arch-clang-15.0_debug named IGSKCI164LM006 by baagaard Wed Jun  5 15:18:41 2024
[0]PETSC ERROR: Configure options --PETSC_ARCH=arch-clang-15.0_debug --with-debugging=1 --with-clanguage=c --with-mpi-compilers=1 --with-shared-libraries=1 --with-64-bit-points=1 --with-large-file-io=1 --with-lgrind=0 --download-chaco=1 --download-parmetis=1 --download-metis=1 --download-triangle --download-ml=1 --download-superlu=1 --with-fc=0 --download-f2cblaslapack --with-hdf5=1 --with-hdf5-include=/software/unix/hdf5-1.14/clang-15.0/include --with-hdf5-lib=/software/unix/hdf5-1.14/clang-15.0/lib/libhdf5.dylib --with-zlib=1
[0]PETSC ERROR: #1 static PetscErrorCode pylith::topology::FieldQuery::queryDBPointFn(PylithInt, PylithReal, const PylithReal *, PylithInt, PylithScalar *, void *)() at /src/cig/pylith/libsrc/pylith/topology/FieldQuery.cc:316
[0]PETSC ERROR: #2 DMProjectPoint_Func_Private() at /software/unix/petsc-dev/src/dm/impls/plex/plexproject.c:128
[0]PETSC ERROR: #3 DMProjectPoint_Private() at /software/unix/petsc-dev/src/dm/impls/plex/plexproject.c:483
[0]PETSC ERROR: #4 DMProjectLocal_Generic_Plex() at /software/unix/petsc-dev/src/dm/impls/plex/plexproject.c:1003
[0]PETSC ERROR: #5 DMProjectFunctionLocal_Plex() at /software/unix/petsc-dev/src/dm/impls/plex/plexproject.c:1034
[0]PETSC ERROR: #6 DMProjectFunctionLocal() at /software/unix/petsc-dev/src/dm/interface/dm.c:8183
[0]PETSC ERROR: #7 void pylith::topology::FieldQuery::queryDB()() at /src/cig/pylith/libsrc/pylith/topology/FieldQuery.cc:223
Fatal error. Calling MPI_Abort() to abort PyLith application.
Traceback (most recent call last):
  File "/software/unix/py3.12-venv/pylith-debug/lib/python3.12/site-packages/pylith/apps/PetscApplication.py", line 55, in onComputeNodes
    self.main(*args, **kwds)
  File "/software/unix/py3.12-venv/pylith-debug/lib/python3.12/site-packages/pylith/apps/PyLithApp.py", line 103, in main
    self.problem.initialize()
  File "/software/unix/py3.12-venv/pylith-debug/lib/python3.12/site-packages/pylith/problems/Problem.py", line 187, in initialize
    ModuleProblem.initialize(self)
  File "/software/unix/py3.12-venv/pylith-debug/lib/python3.12/site-packages/pylith/problems/problems.py", line 165, in initialize
    return _problems.Problem_initialize(self)
           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
RuntimeError: Error detected while in PETSc function.
Abort(-1) on node 0 (rank 0 in comm 0): application called MPI_Abort(MPI_COMM_WORLD, -1) - process 0
/software/unix/py3.12-venv/pylith-debug/bin/nemesis: mpiexec: exit 255
```

## Troubleshooting Strategy

We still get an error about not being able to find an initiation time for a point.
This suggests there are still one or more errors in our spatial database file related to interpolation.
We examine the header and data points for errors.
We notice that our deepest point has a y coordinate of -25 km, but PyLith is looking for values at a point with a y coordinate of -27.621 km.
We need to add an additional point to our spatial database.
This explains why we had `num-locs=4` when we started!

## Resolution

```{code-block} cfg
---
caption: Correct error in `fault_slip.spatialdb`.
---
num-locs = 4
...
0.0   99.0     -2.0       0.0   0.0
0.0  -20.0     -2.0       0.0   0.0
0.0  -25.0      0.0       0.0   0.0
0.0  -99.0      0.0       0.0   0.0
```

Our simulation now runs without errors and the output looks correct.
