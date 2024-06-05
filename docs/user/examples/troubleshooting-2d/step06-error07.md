# Step 6: Error 7

## Error Message

```{code-block} console
---
caption: Error message 7 when running Step 6.
linenos: True
emphasize-lines: 92
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
 -- Pre-initializing fault 'fault'.
 >> /software/unix/py3.12-venv/pylith-debug/lib/python3.12/site-packages/pylith/faults/FaultCohesiveKin.py:87:preinitialize
 -- faultcohesivekin(info)
 -- Pre-initializing fault 'splay'.
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
 -- Pre-initializing fault 'fault'.
 >> /software/unix/py3.12-venv/pylith-debug/lib/python3.12/site-packages/pylith/faults/FaultCohesiveKin.py:87:preinitialize
 -- faultcohesivekin(info)
 -- Pre-initializing fault 'splay'.
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
[0]PETSC ERROR: --------------------- Error Message --------------------------------------------------------------
[0]PETSC ERROR: Residual norm computed by GMRES recursion formula 2.28464e+17 is far from the computed residual norm 2.99942e+27 at restart, residual norm at start of cycle 9.77239e+24
[0]PETSC ERROR: WARNING! There are unused option(s) set! Could be the program crashed before usage or a spelling mistake, etc!
[0]PETSC ERROR:   Option left: name:-ksp_converged_reason (no value) source: code
[0]PETSC ERROR:   Option left: name:-snes_converged_reason (no value) source: code
[0]PETSC ERROR: See https://petsc.org/release/faq/ for trouble shooting.
[0]PETSC ERROR: Petsc Development GIT revision: v3.21.2-167-g4fed2113cae  GIT Date: 2024-05-31 10:11:14 -0400
[0]PETSC ERROR: /software/unix/py3.12-venv/pylith-debug/bin/mpinemesis on a arch-clang-15.0_debug named IGSKCI164LM006 by baagaard Wed Jun  5 13:56:54 2024
[0]PETSC ERROR: Configure options --PETSC_ARCH=arch-clang-15.0_debug --with-debugging=1 --with-clanguage=c --with-mpi-compilers=1 --with-shared-libraries=1 --with-64-bit-points=1 --with-large-file-io=1 --with-lgrind=0 --download-chaco=1 --download-parmetis=1 --download-metis=1 --download-triangle --download-ml=1 --download-superlu=1 --with-fc=0 --download-f2cblaslapack --with-hdf5=1 --with-hdf5-include=/software/unix/hdf5-1.14/clang-15.0/include --with-hdf5-lib=/software/unix/hdf5-1.14/clang-15.0/lib/libhdf5.dylib --with-zlib=1
[0]PETSC ERROR: #1 KSPGMRESCycle() at /software/unix/petsc-dev/src/ksp/ksp/impls/gmres/gmres.c:116
[0]PETSC ERROR: #2 KSPSolve_GMRES() at /software/unix/petsc-dev/src/ksp/ksp/impls/gmres/gmres.c:227
[0]PETSC ERROR: #3 KSPSolve_Private() at /software/unix/petsc-dev/src/ksp/ksp/interface/itfunc.c:905
[0]PETSC ERROR: #4 KSPSolve() at /software/unix/petsc-dev/src/ksp/ksp/interface/itfunc.c:1078
[0]PETSC ERROR: #5 SNESSolve_NEWTONLS() at /software/unix/petsc-dev/src/snes/impls/ls/ls.c:220
[0]PETSC ERROR: #6 SNESSolve() at /software/unix/petsc-dev/src/snes/interface/snes.c:4738
[0]PETSC ERROR: #7 TSTheta_SNESSolve() at /software/unix/petsc-dev/src/ts/impls/implicit/theta/theta.c:174
[0]PETSC ERROR: #8 TSStep_Theta() at /software/unix/petsc-dev/src/ts/impls/implicit/theta/theta.c:225
[0]PETSC ERROR: #9 TSStep() at /software/unix/petsc-dev/src/ts/interface/ts.c:3391
[0]PETSC ERROR: #10 TSSolve() at /software/unix/petsc-dev/src/ts/interface/ts.c:4037
[0]PETSC ERROR: #11 void pylith::problems::TimeDependent::solve()() at /src/cig/pylith/libsrc/pylith/problems/TimeDependent.cc:487
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
/software/unix/py3.12-venv/pylith-debug/bin/pylith: /software/unix/py3.12-venv/pylith-debug/bin/nemesis: exit 1
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
