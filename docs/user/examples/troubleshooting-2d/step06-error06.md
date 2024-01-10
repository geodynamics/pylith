# Step 6: Error 6

## Error Message

```{code-block} console
---
caption: Error message 6 when running Step 6.
linenos: True
emphasize-lines: 9-10
---
$ pylith step06_twofaults.cfg

# -- many lines omitted --

 -- Solving problem.
0 TS dt 0.01 time 0.
    0 SNES Function norm 2.001189838638e-02 
[0]PETSC ERROR: --------------------- Error Message --------------------------------------------------------------
[0]PETSC ERROR: Zero pivot in LU factorization: https://petsc.org/release/faq/#zeropivot
[0]PETSC ERROR: Zero pivot row 78 value 1.11022e-16 tolerance 2.22045e-14
# -- lines with PETSc configuration omitted --
[0]PETSC ERROR: #1 MatPivotCheck_none() at /software/baagaard/petsc-dev/include/petsc/private/matimpl.h:802
[0]PETSC ERROR: #2 MatPivotCheck() at /software/baagaard/petsc-dev/include/petsc/private/matimpl.h:821
# -- many lines of PETSc stack omitted --
[0]PETSC ERROR: #34 void pylith::problems::TimeDependent::deallocate()() at /home/pylith-user/src/cig/pylith/libsrc/pylith/problems/TimeDependent.cc:92
Traceback (most recent call last):
  File "/software/baagaard/py38-venv/pylith-debug/lib/python3.8/site-packages/pylith/apps/PetscApplication.py", line 61, in onComputeNodes
    self.main(*args, **kwds)
  File "/software/baagaard/py38-venv/pylith-debug/lib/python3.8/site-packages/pylith/apps/PyLithApp.py", line 120, in main
    self.problem.run(self)
  File "/software/baagaard/py38-venv/pylith-debug/lib/python3.8/site-packages/pylith/problems/TimeDependent.py", line 141, in run
    ModuleTimeDependent.solve(self)
  File "/software/baagaard/py38-venv/pylith-debug/lib/python3.8/site-packages/pylith/problems/problems.py", line 223, in solve
    return _problems.TimeDependent_solve(self)
RuntimeError: Error detected while in PETSc function.
```

## Troubleshooting Strategy

PETSc encounters a zero pivot during the LU factorization.
This suggests there is an error in how we setup the problem.
We examine the parameters in `step06_twofaults.cfg`, and we do not find any obvious errors.
We load the `info` files into ParaView to visualize their contents for errors.
After loading `step06_twofaults-fault_info.xmf` and `step06_twofaults-splay_info.xmf`, we see that the faults cross each other as shown in {numref}`fig:example:troubleshooting:2d:step06:faults1`.
In `step06_twofaults.cfg` we see that the splay fault is listed first in the array of faults.
The through-going fault (main fault) should be listed first.

:::{figure-md} fig:example:troubleshooting:2d:step06:faults1
<img src="figs/step06-faults-wrong1.*" alt="" scale="75%">

Incorrect geometry for the splay and main faults.
The splay fault crosses the main fault instead of terminating where they intersect.
:::

## Resolution

```{code-block} cfg
---
caption: Correct error in `step06_twofaults.cfg`.
---
# Error
[pylithapp.problem]
interfaces = [splay, fault]

# Correct
[pylithapp.problem]
interfaces = [fault, splay]
```
