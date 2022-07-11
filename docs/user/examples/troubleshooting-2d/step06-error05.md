# Step 6: Error 5

## Error Message

```{code-block} console
---
caption: Error message 5 when running Step 6.
linenos: True
emphasize-lines: 14-17
---
$ pylith step06_twofaults.cfg

 -- Initializing timedependent problem with quasistatic formulation.
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
RuntimeError: Could not find value 'final_slip_opening' in spatial database 'Fault rupture for main fault'. Available values are:
  final-slip-left-lateral
  final-slip-opening
  initiation-time
```

## Troubleshooting Strategy

We have more errors with `fault_slip.spatialdb`.
The error message on lines 14-17 shows that PyLith is looking for `final_slip_opening` in the spatial database, but it found `final-slip-opening` instead.
We need to change the dashes (used in PyLith v1.x and v2.x) to underscores (used in PyLith v3.x); we made this change to be consistent with the names of the output fields.

## Resolution

```{code-block} cfg
---
caption: Correct error in `fault_slip.spatialdb`.
---
# Error
value-names = final-slip-left-lateral  final-slip-opening  initiation-time

# Correct
value-names = final_slip_left_lateral  final_slip_opening  initiation_time
```
