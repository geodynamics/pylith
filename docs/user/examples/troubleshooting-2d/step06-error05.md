# Step 6: Error 5

## Error Message

```{code-block} console
---
caption: Error message 5 when running Step 6.
linenos: True
emphasize-lines: 97-100
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
 >> /software/unix/py3.12-venv/pylith-debug/lib/python3.12/site-packages/pylith/problems/Problem.py:219:_printInfo
 -- timedependent(info)
 -- Scales for nondimensionalization:
    Length scale: 1000*m
    Time scale: 3.15576e+09*s
    Pressure scale: 3e+10*m**-1*kg*s**-2
    Density scale: 2.98765e+23*m**-3*kg
    Temperature scale: 1*K
 >> /software/unix/py3.12-venv/pylith-debug/lib/python3.12/site-packages/pylith/problems/Problem.py:185:initialize
 -- timedependent(info)
 -- Initializing timedependent problem with quasistatic formulation.
 >> /src/cig/pylith/libsrc/pylith/utils/PetscOptions.cc:239:static void pylith::utils::_PetscOptions::write(pythia::journal::info_t &, const char *, const PetscOptions &)
 -- petscoptions(info)
 -- Setting PETSc options:
dm_reorder_section = true
dm_reorder_section_type = cohesive
ksp_atol = 1.0e-12
ksp_converged_reason = true
ksp_error_if_not_converged = true
ksp_guess_pod_size = 8
ksp_guess_type = pod
ksp_rtol = 1.0e-12
mg_fine_pc_type = vpbjacobi
pc_type = gamg
snes_atol = 1.0e-9
snes_converged_reason = true
snes_error_if_not_converged = true
snes_monitor = true
snes_rtol = 1.0e-12
ts_error_if_step_fails = true
ts_monitor = true
ts_type = beuler

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
RuntimeError: Could not find value 'final_slip_opening' in spatial database 'Fault rupture for main fault'. Available values are:
  final-slip-left-lateral
  final-slip-opening
  initiation-time

/software/unix/py3.12-venv/pylith-debug/bin/nemesis: mpiexec: exit 255
/software/unix/py3.12-venv/pylith-debug/bin/pylith: /software/unix/py3.12-venv/pylith-debug/bin/nemesis: exit 1
```

## Troubleshooting Strategy

We have more errors with `fault_slip.spatialdb`.
The error message on lines 97-100 shows that PyLith is looking for `final_slip_opening` in the spatial database, but it found `final-slip-opening` instead.
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
