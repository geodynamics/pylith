# Step 6: Error 7

## Error Message

```{code-block} console
---
caption: Error message 7 when running Step 6.
linenos: True
emphasize-lines: 9-10
---
$ pylith step06_twofaults.cfg

# -- many lines omitted --

 -- Solving problem.
0 TS dt 0.01 time 0.
    0 SNES Function norm 2.024143393875e-02 
[0]PETSC ERROR: --------------------- Error Message --------------------------------------------------------------
[0]PETSC ERROR: Residual norm computed by GMRES recursion formula 3.48613e+10 is far from the computed residual norm 6.92443e+12 at restart, residual norm at start of cycle 6.91369e+12
# -- lines with PETSc configuration omitted --
[0]PETSC ERROR: #1 KSPGMRESCycle() at /software/baagaard/petsc-dev/src/ksp/ksp/impls/gmres/gmres.c:126
[0]PETSC ERROR: #2 KSPSolve_GMRES() at /software/baagaard/petsc-dev/src/ksp/ksp/impls/gmres/gmres.c:243
# -- many lines of PETSc stack omitted --
[0]PETSC ERROR: #11 void pylith::problems::TimeDependent::solve()() at /home/pylith-user/src/cig/pylith/libsrc/pylith/problems/TimeDependent.cc:429
Fatal error. Calling MPI_Abort() to abort PyLith application.
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

PETSc encounters another error during the solve.
This also suggests there is an error in how we setup the problem.
We again load the `info` files into ParaView to visualize their contents for errors.
After loading `step06_twofaults-fault_info.xmf` and `step06_twofaults-splay_info.xmf`, we see that the faults overlap just below their intersection as shown in {numref}`fig:example:troubleshooting:2d:step06:faults2`.
The bottom of each fault is buried, so we need to identify the buried edges so that when PyLith inserts the cohesive cells it can terminate the fault properly.

:::{figure-md} fig:example:troubleshooting:2d:step06:faults2
<img src="figs/step06-faults-wrong2.*" alt="" scale="75%">

Incorrect geometry for the splay and main faults.
The two faults overlap just below their intersection.
The bottom edge of each fault is buried, so we must identify the buried edges so that PyLith properly terminates the edges of the faults when inserting the cohesive cells.
:::

## Resolution

```{code-block} cfg
---
caption: Correct error in `step06_twofaults.cfg`.
---
[pylithapp.problem.interfaces.fault]
...
edge = fault_end
edge_value = 21

[pylithapp.problem.interfaces.splay]
...
edge = splay_end
edge_value = 23
```
