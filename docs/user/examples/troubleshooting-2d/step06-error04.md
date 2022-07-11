# Step 6: Error 4

## Error Message

```{code-block} console
---
caption: Error message 4 when running Step 6.
linenos: True
emphasize-lines: 14-15
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
RuntimeError: Error occurred while reading spatial database file 'fault_slip.spatialdb'.
Read data for 3 out of 4 points.
Error reading coordinates from buffer ''.
```

## Troubleshooting Strategy

The error message on lines 14-15 indicates there is an error reading the `fault_slip.spatialdb` spatial database for the fault slip.
PyLith was able to read data for 3 of 4 points.
The file `fault_slip.spatialdb` contains only 3 points but `num-locs` is 4.

## Resolution

```{code-block} cfg
---
caption: Correct error in `fault_slip.spatialdb`.
---
# Error
num-locs = 4

# Correct
num-locs = 3
```
