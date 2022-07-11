# Step 1: Error 4

## Error Message

```{code-block} console
---
caption: Error message 4 when running Step 1.
linenos: True
emphasize-lines: 20
---
$ pylith step01_gravity.cfg

 >> /software/baagaard/py38-venv/pylith-debug/lib/python3.8/site-packages/pylith/meshio/MeshIOObj.py:44:read
 -- meshiopetsc(info)
 -- Reading finite-element mesh

# -- many lines omitted --

 -- Verifying compatibility of problem configuration.
Fatal error. Calling MPI_Abort() to abort PyLith application.
Traceback (most recent call last):
  File "/software/baagaard/py38-venv/pylith-debug/lib/python3.8/site-packages/pylith/apps/PetscApplication.py", line 61, in onComputeNodes
    self.main(*args, **kwds)
  File "/software/baagaard/py38-venv/pylith-debug/lib/python3.8/site-packages/pylith/apps/PyLithApp.py", line 108, in main
    self.problem.verifyConfiguration()
  File "/software/baagaard/py38-venv/pylith-debug/lib/python3.8/site-packages/pylith/problems/Problem.py", line 177, in verifyConfiguration
    ModuleProblem.verifyConfiguration(self)
  File "/software/baagaard/py38-venv/pylith-debug/lib/python3.8/site-packages/pylith/problems/problems.py", line 167, in verifyConfiguration
    return _problems.Problem_verifyConfiguration(self)
RuntimeError: Material label_value '3' for cell '3009' does not match the label_value of any materials or interfaces.
```

## Troubleshooting Strategy

The output shows a Python Traceback and then the error message on Line 20.
The error indicates the finite-element mesh file contains a cell with a label value of 3, but the parameter files do not have a material with a label value of 3.
We examine the `pylithapp.problem.materials` sections of `pylithapp.cfg` and see that the label values are 0, 1, and 2 in the parameter file rather than 1, 2, and 3.

## Resolution

```{code-block} cfg
---
caption: Correct error in `pylithapp.cfg`.
---
[pylithapp.problem.materials.slab]
label_value = 1
...

[pylithapp.problem.materials.crust]
label_value = 2
...

[pylithapp.problem.materials.wedge]
label_value = 3
...
```
