# Step 6: Error 9

## Error Message

```{code-block} console
---
caption: Error message 9 when running Step 6.
linenos: True
emphasize-lines: 6
---
$ pylith step06_twofaults.cfg

 -- Initializing timedependent problem with quasistatic formulation.
[0]PETSC ERROR: --------------------- Error Message --------------------------------------------------------------
[0]PETSC ERROR: Error in external library
[0]PETSC ERROR: Could not find values for initiation_time at (  -24329  -29046.3) in spatial database 'Fault rupture for main fault'.
# -- lines with PETSc configuration omitted --
[0]PETSC ERROR: #1 static PetscErrorCode pylith::topology::FieldQuery::queryDBPointFn(PylithInt, PylithReal, const PylithReal*, PylithInt, PylithScalar*, void*)() at /home/baagaard/src/cig/pylith/libsrc/pylith/topology/FieldQuery.cc:313
# -- many lines of PETSc stack omitted --
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

With the linear interpolation we get an error about not being able to find an initiation time for a point.
This suggests there are is one or more errors in our spatial database file related to interpolation.
We examine the header and data points for errors.
We notice that our points lie along a line (data dimension is 1), but our header has `data-dim=2`.

## Resolution

```{code-block} cfg
---
caption: Correct error in `step06_twofaults.cfg`.
---
# Error
data-dim = 2
    
# Correct
data-dim = 1
```

Our simulation now runs without errors and the output looks correct.
