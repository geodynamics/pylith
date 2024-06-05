# Step 6: Error 3

## Error Message

```{code-block} console
---
caption: Error message 3 when running Step 6.
linenos: True
emphasize-lines: 64-65
---
$ pylith step06_twofaults.cfg

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
 >> /software/unix/py3.12-venv/pylith-debug/lib/python3.12/site-packages/pylith/faults/FaultCohesiveKin.py:87:preinitialize
 -- faultcohesivekin(info)
 -- Pre-initializing fault 'splay'.
 >> /software/unix/py3.12-venv/pylith-debug/lib/python3.12/site-packages/pylith/faults/FaultCohesiveKin.py:87:preinitialize
 -- faultcohesivekin(info)
 -- Pre-initializing fault 'fault'.
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
 >> /software/unix/py3.12-venv/pylith-debug/lib/python3.12/site-packages/pylith/faults/FaultCohesiveKin.py:87:preinitialize
 -- faultcohesivekin(info)
 -- Pre-initializing fault 'splay'.
 >> /software/unix/py3.12-venv/pylith-debug/lib/python3.12/site-packages/pylith/faults/FaultCohesiveKin.py:87:preinitialize
 -- faultcohesivekin(info)
 -- Pre-initializing fault 'fault'.
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
RuntimeError: Could not find 'lagrange_multiplier_fault' in solution field. Field contains: 'displacement'Abort(-1) on node 0 (rank 0 in comm 0): application called MPI_Abort(MPI_COMM_WORLD, -1) - process 0
; the missing fields are required for interface 'FaultCohesiveImpulses'.
/software/unix/py3.12-venv/pylith-debug/bin/nemesis: mpiexec: exit 255
/software/unix/py3.12-venv/pylith-debug/bin/pylith: /software/unix/py3.12-venv/pylith-debug/bin/nemesis: exit 1
```

## Troubleshooting Strategy

We no longer have errors during the problem configuration.
Now we have errors while doing additional verification of the problem.
After the Python Traceback, we see the error message on lines 64-65.
The faults check to make sure the solution field contains the necessary subfields, and the splay fault cannot find the `lagrange_multiplier_fault` subfield.
The easiest way to diagnose an error like this is to view the JSON file automatically generated by PyLith; it contains all of the parameters, including any defaults used.
We point our web browser to <https://geodynamics.github.io/pylith_parameters/> and load the parameter file `output/step06_twofaults-parameters.json`.
In the left panel we navigate to the solution field and see that it the subfields are set to `pylith.problems.SolnDisp`, so that the solution field only contains a single subfield, `displacement`.
We want the solution field to contain both `displacement` and `lagrange_multiplier_fault`.

## Resolution

```{code-block} cfg
---
caption: Correct error in `step06_twofaults.cfg`.
---
[pylithapp.problem]
solution = pylith.problems.SolnDispLagrange
```
