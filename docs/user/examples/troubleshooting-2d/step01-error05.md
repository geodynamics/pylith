# Step 1: Error 5

## Error Message

```{code-block} console
---
caption: Error message 5 when running Step 1.
linenos: True
emphasize-lines: 8,9
---
$ pylith step01_gravity.cfg

 >> /software/baagaard/py38-venv/pylith-debug/lib/python3.8/site-packages/pylith/problems/Problem.py:186:initialize
 -- timedependent(info)
 -- Initializing timedependent problem with quasistatic formulation.
[0]PETSC ERROR: --------------------- Error Message --------------------------------------------------------------
[0]PETSC ERROR: Error in external library
[0]PETSC ERROR: Error converting spatial database values for gravitational_acceleration at (  -96623.5  -72650.4) in
spatial database 'Gravity field'. Found near zero magnitude (0) for gravity field vector (  0  0).

# -- lines with PETSc configuration omitted --

[0]PETSC ERROR: #1 static PetscErrorCode pylith::topology::FieldQuery::queryDBPointFn(PylithInt, PylithReal, const PylithReal*, PylithInt, PylithScalar*, void*)() at /home/baagaard/src/cig/pylith/libsrc/pylith/topology/FieldQuery.cc:327
[0]PETSC ERROR: #2 DMProjectPoint_Func_Private() at /software/baagaard/petsc-dev/src/dm/impls/plex/plexproject.c:127
[0]PETSC ERROR: #3 DMProjectPoint_Private() at /software/baagaard/petsc-dev/src/dm/impls/plex/plexproject.c:409
[0]PETSC ERROR: #4 DMProjectLocal_Generic_Plex() at /software/baagaard/petsc-dev/src/dm/impls/plex/plexproject.c:902
[0]PETSC ERROR: #5 DMProjectFunctionLocal_Plex() at /software/baagaard/petsc-dev/src/dm/impls/plex/plexproject.c:933
[0]PETSC ERROR: #6 DMProjectFunctionLocal() at /software/baagaard/petsc-dev/src/dm/interface/dm.c:8869
[0]PETSC ERROR: #7 void pylith::topology::FieldQuery::queryDB()() at /home/baagaard/src/cig/pylith/libsrc/pylith/topology/FieldQuery.cc:211
Fatal error. Calling MPI_Abort() to abort PyLith application.
Traceback (most recent call last):
  File "/software/baagaard/py38-venv/pylith-debug/lib/python3.8/site-packages/pylith/apps/PetscApplication.py", line 61, in onComputeNodes
    self.main(*args, **kwds)
  File "/software/baagaard/py38-venv/pylith-debug/lib/python3.8/site-packages/pylith/apps/PyLithApp.py", line 110, in main
    self.problem.initialize()
  File "/software/baagaard/py38-venv/pylith-debug/lib/python3.8/site-packages/pylith/problems/Problem.py", line 188, in initialize
    ModuleProblem.initialize(self)
  File "/software/baagaard/py38-venv/pylith-debug/lib/python3.8/site-packages/pylith/problems/problems.py", line 170, in initialize
    return _problems.Problem_initialize(self)
RuntimeError: Error detected while in PETSc function.
```

## Troubleshooting Strategy

During initialization of the time-dependent problem, PyLith encountered an error while calling PETSc functions.
When an error occurs within a PETSc routine, PETSc dumps a lot of information to the screen.
First, it display the error message; second, it displays configuration information; third, it displays the stacktrace.
PyLith traps the PETSc error and displays the Python Traceback.
In this case, the error message on lines 8-9 indicate that the gravity field vector has a magnitude of zero.
We notice that it shows only the first two components because this is a 2D problem, whereas the [default](https://spatialdata.readthedocs.io/en/latest/user/components/spatialdb/GravityField.html) is (0, 0, -1) and intended for 3D problems.

## Resolution

```{code-block} cfg
---
caption: Correct error in `step01_gravity.cfg`.
---
[pylithapp.problem]
gravity_field = spatialdata.spatialdb.GravityField
gravity_field.gravity_dir = [0.0, -1.0, 0.0]
```

Our simulation now runs without errors and the output looks correct.
