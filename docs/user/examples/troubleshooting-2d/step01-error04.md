# Step 1: Error 4

## Error Message

```{code-block} console
---
caption: Error message 4 when running Step 1.
linenos: True
emphasize-lines: 52
---
$ pylith step01a_gravity.cfg

 >> /software/unix/py3.12-venv/pylith-debug/lib/python3.12/site-packages/pylith/apps/PyLithApp.py:77:main
 -- pylithapp(info)
 -- Running on 1 process(es).
 >> /software/unix/py3.12-venv/pylith-debug/lib/python3.12/site-packages/pylith/meshio/MeshIOObj.py:38:read
 -- meshiopetsc(info)
 -- Reading finite-element mesh
 >> /src/cig/pylith/libsrc/pylith/meshio/MeshIO.cc:85:void pylith::meshio::MeshIO::read(pylith::topology::Mesh *, const bool)
 -- meshiopetsc(info)
 -- Component 'reader': Domain bounding box:
    (-100000, 100000)
    (-100000, 0)
 >> /software/unix/py3.12-venv/pylith-debug/lib/python3.12/site-packages/pylith/problems/Problem.py:116:preinitialize
 -- timedependent(info)
 -- Performing minimal initialization before verifying configuration.
 >> /software/unix/py3.12-venv/pylith-debug/lib/python3.12/site-packages/pylith/problems/Solution.py:39:preinitialize
 -- solution(info)
 -- Performing minimal initialization of solution.
 >> /software/unix/py3.12-venv/pylith-debug/lib/python3.12/site-packages/pylith/materials/RheologyElasticity.py:35:preinitialize
 -- isotropiclinearelasticity(info)
 -- Performing minimal initialization of elasticity rheology 'bulk_rheology'.
 >> /software/unix/py3.12-venv/pylith-debug/lib/python3.12/site-packages/pylith/materials/RheologyElasticity.py:35:preinitialize
 -- isotropiclinearelasticity(info)
 -- Performing minimal initialization of elasticity rheology 'bulk_rheology'.
 >> /software/unix/py3.12-venv/pylith-debug/lib/python3.12/site-packages/pylith/materials/RheologyElasticity.py:35:preinitialize
 -- isotropiclinearelasticity(info)
 -- Performing minimal initialization of elasticity rheology 'bulk_rheology'.
 >> /software/unix/py3.12-venv/pylith-debug/lib/python3.12/site-packages/pylith/bc/DirichletTimeDependent.py:86:preinitialize
 -- dirichlettimedependent(info)
 -- Performing minimal initialization of time-dependent Dirichlet boundary condition 'bc_xneg'.
 >> /software/unix/py3.12-venv/pylith-debug/lib/python3.12/site-packages/pylith/bc/DirichletTimeDependent.py:86:preinitialize
 -- dirichlettimedependent(info)
 -- Performing minimal initialization of time-dependent Dirichlet boundary condition 'bc_xpos'.
 >> /software/unix/py3.12-venv/pylith-debug/lib/python3.12/site-packages/pylith/bc/DirichletTimeDependent.py:86:preinitialize
 -- dirichlettimedependent(info)
 -- Performing minimal initialization of time-dependent Dirichlet boundary condition 'bc_yneg'.
 >> /software/unix/py3.12-venv/pylith-debug/lib/python3.12/site-packages/pylith/problems/Problem.py:174:verifyConfiguration
 -- timedependent(info)
 -- Verifying compatibility of problem configuration.
Fatal error. Calling MPI_Abort() to abort PyLith application.
Traceback (most recent call last):
  File "/software/unix/py3.12-venv/pylith-debug/lib/python3.12/site-packages/pylith/apps/PetscApplication.py", line 55, in onComputeNodes
    self.main(*args, **kwds)
  File "/software/unix/py3.12-venv/pylith-debug/lib/python3.12/site-packages/pylith/apps/PyLithApp.py", line 101, in main
    self.problem.verifyConfiguration()
  File "/software/unix/py3.12-venv/pylith-debug/lib/python3.12/site-packages/pylith/problems/Problem.py", line 176, in verifyConfiguration
    ModuleProblem.verifyConfiguration(self)
  File "/software/unix/py3.12-venv/pylith-debug/lib/python3.12/site-packages/pylith/problems/problems.py", line 162, in verifyConfiguration
    return _problems.Problem_verifyConfiguration(self)
           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
RuntimeError: Material label_value '3' for cell '609' does not match the label_value of any materials or interfaces.
Abort(-1) on node 0 (rank 0 in comm 0): application called MPI_Abort(MPI_COMM_WORLD, -1) - process 0
/software/unix/py3.12-venv/pylith-debug/bin/nemesis: mpiexec: exit 255
/software/unix/py3.12-venv/pylith-debug/bin/pylith: /software/unix/py3.12-venv/pylith-debug/bin/nemesis: exit 1
```

## Troubleshooting Strategy

The output shows a Python Traceback and then the error message on Line 52.
The error indicates the finite-element mesh file contains a cell with a label value of 3, but the parameter files do not have a material with a label value of 3.
We examine the `pylithapp.problem.materials` sections of `pylithapp.cfg` and see that the label values are 0, 1, and 2 in the parameter file rather than 1, 2, and 3 that we set in the Gmsh Python script.

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
