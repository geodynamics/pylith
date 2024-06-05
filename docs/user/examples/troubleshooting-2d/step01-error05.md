# Step 1: Error 5

## Error Message

```{code-block} console
---
caption: Error message 5 when running Step 1.
linenos: True
emphasize-lines: 73-74
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
ksp_atol = 1.0e-12
ksp_converged_reason = true
ksp_error_if_not_converged = true
ksp_guess_pod_size = 8
ksp_guess_type = pod
ksp_rtol = 1.0e-12
pc_type = lu
snes_atol = 1.0e-9
snes_converged_reason = true
snes_error_if_not_converged = true
snes_monitor = true
snes_rtol = 1.0e-12
ts_error_if_step_fails = true
ts_monitor = true
ts_type = beuler

[0]PETSC ERROR: --------------------- Error Message --------------------------------------------------------------
[0]PETSC ERROR: Error in external library
[0]PETSC ERROR: Error converting spatial database values for gravitational_acceleration at (  -92285.3  -72772.7) in spatial database 'Gravity field'. Found near zero magnitude (0) for gravity field vector (  0  0).
[0]PETSC ERROR: WARNING! There are unused option(s) set! Could be the program crashed before usage or a spelling mistake, etc!
[0]PETSC ERROR:   Option left: name:-ksp_atol value: 1.0e-12 source: code
[0]PETSC ERROR:   Option left: name:-ksp_converged_reason (no value) source: code
[0]PETSC ERROR:   Option left: name:-ksp_error_if_not_converged (no value) source: code
[0]PETSC ERROR:   Option left: name:-ksp_guess_pod_size value: 8 source: code
[0]PETSC ERROR:   Option left: name:-ksp_guess_type value: pod source: code
[0]PETSC ERROR:   Option left: name:-ksp_rtol value: 1.0e-12 source: code
[0]PETSC ERROR:   Option left: name:-pc_type value: lu source: code
[0]PETSC ERROR:   Option left: name:-snes_atol value: 1.0e-9 source: code
[0]PETSC ERROR:   Option left: name:-snes_converged_reason (no value) source: code
[0]PETSC ERROR:   Option left: name:-snes_error_if_not_converged (no value) source: code
[0]PETSC ERROR:   Option left: name:-snes_monitor (no value) source: code
[0]PETSC ERROR:   Option left: name:-snes_rtol value: 1.0e-12 source: code
[0]PETSC ERROR:   Option left: name:-ts_error_if_step_fails (no value) source: code
[0]PETSC ERROR:   Option left: name:-ts_monitor (no value) source: code
[0]PETSC ERROR:   Option left: name:-ts_type value: beuler source: code
[0]PETSC ERROR: See https://petsc.org/release/faq/ for trouble shooting.
[0]PETSC ERROR: Petsc Development GIT revision: v3.21.2-167-g4fed2113cae  GIT Date: 2024-05-31 10:11:14 -0400
[0]PETSC ERROR: /software/unix/py3.12-venv/pylith-debug/bin/mpinemesis on a arch-clang-15.0_debug named IGSKCI164LM006 by baagaard Wed Jun  5 13:18:39 2024
[0]PETSC ERROR: Configure options --PETSC_ARCH=arch-clang-15.0_debug --with-debugging=1 --with-clanguage=c --with-mpi-compilers=1 --with-shared-libraries=1 --with-64-bit-points=1 --with-large-file-io=1 --with-lgrind=0 --download-chaco=1 --download-parmetis=1 --download-metis=1 --download-triangle --download-ml=1 --download-superlu=1 --with-fc=0 --download-f2cblaslapack --with-hdf5=1 --with-hdf5-include=/software/unix/hdf5-1.14/clang-15.0/include --with-hdf5-lib=/software/unix/hdf5-1.14/clang-15.0/lib/libhdf5.dylib --with-zlib=1
[0]PETSC ERROR: #1 static PetscErrorCode pylith::topology::FieldQuery::queryDBPointFn(PylithInt, PylithReal, const PylithReal *, PylithInt, PylithScalar *, void *)() at /src/cig/pylith/libsrc/pylith/topology/FieldQuery.cc:330
[0]PETSC ERROR: #2 DMProjectPoint_Func_Private() at /software/unix/petsc-dev/src/dm/impls/plex/plexproject.c:128
[0]PETSC ERROR: #3 DMProjectPoint_Private() at /software/unix/petsc-dev/src/dm/impls/plex/plexproject.c:483
[0]PETSC ERROR: #4 DMProjectLocal_Generic_Plex() at /software/unix/petsc-dev/src/dm/impls/plex/plexproject.c:1003
[0]PETSC ERROR: #5 DMProjectFunctionLocal_Plex() at /software/unix/petsc-dev/src/dm/impls/plex/plexproject.c:1034
[0]PETSC ERROR: #6 DMProjectFunctionLocal() at /software/unix/petsc-dev/src/dm/interface/dm.c:8183
[0]PETSC ERROR: #7 void pylith::topology::FieldQuery::queryDB()() at /src/cig/pylith/libsrc/pylith/topology/FieldQuery.cc:223
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
RuntimeError: Error detected while in PETSc function.
Abort(-1) on node 0 (rank 0 in comm 0): application called MPI_Abort(MPI_COMM_WORLD, -1) - process 0
/software/unix/py3.12-venv/pylith-debug/bin/nemesis: mpiexec: exit 255
/software/unix/py3.12-venv/pylith-debug/bin/pylith: /software/unix/py3.12-venv/pylith-debug/bin/nemesis: exit 1

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
