# Step 6: Error 6

## Error Message

```{code-block} console
---
caption: Error message 6 when running Step 6.
linenos: True
emphasize-lines: 93
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

 >> /software/unix/py3.12-venv/pylith-debug/lib/python3.12/site-packages/pylith/problems/TimeDependent.py:132:run
 -- timedependent(info)
 -- Solving problem.
0 TS dt 0.2 time -0.2
    0 SNES Function norm 1.989081034695e-02
      Linear solve converged due to CONVERGED_RTOL iterations 2
    1 SNES Function norm 1.989081034695e-02
    Nonlinear solve did not converge due to DIVERGED_MAX_IT iterations 1
[0]PETSC ERROR: --------------------- Error Message --------------------------------------------------------------
[0]PETSC ERROR: SNESSolve has not converged
[0]PETSC ERROR: See https://petsc.org/release/faq/ for trouble shooting.
[0]PETSC ERROR: Petsc Development GIT revision: v3.21.2-167-g4fed2113cae  GIT Date: 2024-05-31 10:11:14 -0400
[0]PETSC ERROR: /software/unix/py3.12-venv/pylith-debug/bin/mpinemesis on a arch-clang-15.0_debug named IGSKCI164LM006 by baagaard Wed Jun  5 13:47:43 2024
[0]PETSC ERROR: Configure options --PETSC_ARCH=arch-clang-15.0_debug --with-debugging=1 --with-clanguage=c --with-mpi-compilers=1 --with-shared-libraries=1 --with-64-bit-points=1 --with-large-file-io=1 --with-lgrind=0 --download-chaco=1 --download-parmetis=1 --download-metis=1 --download-triangle --download-ml=1 --download-superlu=1 --with-fc=0 --download-f2cblaslapack --with-hdf5=1 --with-hdf5-include=/software/unix/hdf5-1.14/clang-15.0/include --with-hdf5-lib=/software/unix/hdf5-1.14/clang-15.0/lib/libhdf5.dylib --with-zlib=1
[0]PETSC ERROR: #1 SNESSolve() at /software/unix/petsc-dev/src/snes/interface/snes.c:4751
[0]PETSC ERROR: #2 TSTheta_SNESSolve() at /software/unix/petsc-dev/src/ts/impls/implicit/theta/theta.c:174
[0]PETSC ERROR: #3 TSStep_Theta() at /software/unix/petsc-dev/src/ts/impls/implicit/theta/theta.c:225
[0]PETSC ERROR: #4 TSStep() at /software/unix/petsc-dev/src/ts/interface/ts.c:3391
[0]PETSC ERROR: #5 TSSolve() at /software/unix/petsc-dev/src/ts/interface/ts.c:4037
[0]PETSC ERROR: #6 void pylith::problems::TimeDependent::solve()() at /src/cig/pylith/libsrc/pylith/problems/TimeDependent.cc:487
Fatal error. Calling MPI_Abort() to abort PyLith application.
Traceback (most recent call last):
  File "/software/unix/py3.12-venv/pylith-debug/lib/python3.12/site-packages/pylith/apps/PetscApplication.py", line 55, in onComputeNodes
    self.main(*args, **kwds)
  File "/software/unix/py3.12-venv/pylith-debug/lib/python3.12/site-packages/pylith/apps/PyLithApp.py", line 114, in main
    self.problem.run(self)
  File "/software/unix/py3.12-venv/pylith-debug/lib/python3.12/site-packages/pylith/problems/TimeDependent.py", line 134, in run
    ModuleTimeDependent.solve(self)
  File "/software/unix/py3.12-venv/pylith-debug/lib/python3.12/site-packages/pylith/problems/problems.py", line 217, in solve
    return _problems.TimeDependent_solve(self)
           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
RuntimeError: Error detected while in PETSc function.
Abort(-1) on node 0 (rank 0 in comm 0): application called MPI_Abort(MPI_COMM_WORLD, -1) - process 0
/software/unix/py3.12-venv/pylith-debug/bin/nemesis: mpiexec: exit 255
```

## Troubleshooting Strategy

The PETSc nonlinear solver fails to converge in 1 iteration.
Because this is a linear problem, we set the maximum number of nonlinear solver (SNES) iterations to 1, so that we would get this error as a way of identifying errors in the problem setup.
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
