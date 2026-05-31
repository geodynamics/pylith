# Step 6: Error 8

## Error Message

```{code-block} console
---
caption: Error message 8 when running Step 6.
linenos: True
emphasize-lines: 55-56
---
$ pylith step06_twofaults.cfg

 >> software/pylith-debug/lib/python3.12/site-packages/pylith/apps/PyLithApp.py:76:main
 -- info (application-flow)
 -- Running on 1 process(es).
 >> src/cig/pylith/libsrc/pylith/utils/PetscOptions.cc:251:static void pylith::utils::_PetscOptions::write(pythia::journal::info_t&, const char*, const pylith::utils::PetscOptions&)
 -- info (application-flow)
 -- Setting PETSc options:
    dm_reorder_section = true
    dm_reorder_section_type = cohesive
    ksp_atol = 1.0e-7
    ksp_converged_reason = true
    ksp_error_if_not_converged = true
    ksp_gmres_restart = 100
    ksp_guess_pod_size = 8
    ksp_guess_type = pod
    ksp_rtol = 1.0e-14
    mg_fine_ksp_max_it = 5
    mg_fine_pc_type = vpbjacobi
    pc_type = gamg
    snes_atol = 5.0e-7
    snes_converged_reason = true
    snes_error_if_not_converged = true
    snes_monitor = true
    snes_rtol = 1.0e-14
    ts_error_if_step_fails = true
    ts_exact_final_time = matchstep
    ts_monitor = true
    ts_type = beuler
    viewer_hdf5_collective = true

 >> src/cig/pylith/libsrc/pylith/meshio/MeshIOPetsc.cc:204:virtual void pylith::meshio::MeshIOPetsc::_read()
 -- info (application-flow)
 -- Component 'meshiopetsc.reader': Reading finite-element mesh from 'mesh_tri.msh'.
 >> src/cig/pylith/libsrc/pylith/meshio/MeshIO.cc:76:void pylith::meshio::MeshIO::read(pylith::topology::Mesh*, bool)
 -- info (application-flow)
 -- Component 'meshiopetsc.reader': Domain bounding box:
    (-100000, 100000)
    (-100000, 0)
 >> src/cig/pylith/libsrc/pylith/problems/TimeDependent.cc:316:virtual void pylith::problems::TimeDependent::verifyConfiguration() const
 -- info (application-flow)
 -- Component 'timedependent.problem': Verifying problem configuration.
 >> software/pylith-debug/lib/python3.12/site-packages/pylith/problems/Problem.py:238:_printInfo
 -- info (application-flow)
 -- Scales for nondimensionalization:
    Length scale: 2500*m
    Displacement scale: 1*m
    Time scale: 3.15576e+09*s
    Rigidity scale: 1e+10*m**-1*kg*s**-2
    Temperature scale: 1*K
 >> src/cig/pylith/libsrc/pylith/problems/TimeDependent.cc:342:virtual void pylith::problems::TimeDependent::initialize()
 -- info (application-flow)
 -- Component 'timedependent.problem': Initializing problem.
[0]PETSC ERROR: --------------------- Error Message --------------------------------------------------------------
[0]PETSC ERROR: Error in external library
[0]PETSC ERROR: Could not find values for initiation_time at (  -21860.3  -27621) in spatial database 'Fault rupture for main fault'.
[0]PETSC ERROR: WARNING! There are unused option(s) set! Could be the program crashed before usage or a spelling mistake, etc!
[0]PETSC ERROR:   Option left: name:-ksp_atol value: 1.0e-7 source: code
[0]PETSC ERROR:   Option left: name:-ksp_converged_reason (no value) source: code
[0]PETSC ERROR:   Option left: name:-ksp_error_if_not_converged (no value) source: code
[0]PETSC ERROR:   Option left: name:-ksp_gmres_restart value: 100 source: code
[0]PETSC ERROR:   Option left: name:-ksp_guess_pod_size value: 8 source: code
[0]PETSC ERROR:   Option left: name:-ksp_guess_type value: pod source: code
[0]PETSC ERROR:   Option left: name:-ksp_rtol value: 1.0e-14 source: code
[0]PETSC ERROR:   Option left: name:-mg_fine_ksp_max_it value: 5 source: code
[0]PETSC ERROR:   Option left: name:-mg_fine_pc_type value: vpbjacobi source: code
[0]PETSC ERROR:   Option left: name:-pc_type value: gamg source: code
[0]PETSC ERROR:   Option left: name:-snes_atol value: 5.0e-7 source: code
[0]PETSC ERROR:   Option left: name:-snes_converged_reason (no value) source: code
[0]PETSC ERROR:   Option left: name:-snes_error_if_not_converged (no value) source: code
[0]PETSC ERROR:   Option left: name:-snes_monitor (no value) source: code
[0]PETSC ERROR:   Option left: name:-snes_rtol value: 1.0e-14 source: code
[0]PETSC ERROR:   Option left: name:-ts_error_if_step_fails (no value) source: code
[0]PETSC ERROR:   Option left: name:-ts_exact_final_time value: matchstep source: code
[0]PETSC ERROR:   Option left: name:-ts_monitor (no value) source: code
[0]PETSC ERROR:   Option left: name:-ts_type value: beuler source: code
[0]PETSC ERROR: See https://petsc.org/release/faq/ for trouble shooting.
[0]PETSC ERROR: PETSc Development Git Revision: v3.25.1-168-g70613f6caab Git Date: 2026-05-21 15:15:28 +0000
[0]PETSC ERROR: software/pylith-debug/bin/mpinemesis with 1 MPI process(es) and PETSC_ARCH arch-pylith-debug on igskci164warlng.gs.doi.net by baagaard Tue May 26 12:55:27 2026
[0]PETSC ERROR: Configure options: --PETSC_ARCH=arch-pylith-debug --with-debugging=1 --with-clanguage=c --with-mpi-compilers=1 --with-shared-libraries=1 --with-64-bit-points=1 --with-large-file-io=1 --with-lgrind=0 --download-parmetis=1 --download-metis=1 --download-triangle --download-ml=1 --download-superlu=1 --with-fc=0 --download-f2cblaslapack --with-hdf5=1 --with-hdf5-dir=software/pylith-debug --with-zlib=1 CFLAGS+=-g
[0]PETSC ERROR: #1 static PetscErrorCode pylith::topology::FieldQuery::queryDBPointFn(PylithInt, PylithReal, const PylithReal*, PylithInt, PylithScalar*, void*)() at src/cig/pylith/libsrc/pylith/topology/FieldQuery.cc:320
[0]PETSC ERROR: #2 DMProjectPoint_Func_Private() at /software/baagaard/petsc-dev/src/dm/impls/plex/plexproject.c:143
[0]PETSC ERROR: #3 DMProjectPoint_Private() at /software/baagaard/petsc-dev/src/dm/impls/plex/plexproject.c:545
[0]PETSC ERROR: #4 DMProjectLocal_Generic_Plex() at /software/baagaard/petsc-dev/src/dm/impls/plex/plexproject.c:1083
[0]PETSC ERROR: #5 DMProjectFunctionLocal_Plex() at /software/baagaard/petsc-dev/src/dm/impls/plex/plexproject.c:1114
[0]PETSC ERROR: #6 DMProjectFunctionLocal() at /software/baagaard/petsc-dev/src/dm/interface/dm.c:8361
[0]PETSC ERROR: #7 void pylith::topology::FieldQuery::queryDB()() at src/cig/pylith/libsrc/pylith/topology/FieldQuery.cc:227
Fatal error. Calling MPI_Abort() to abort PyLith application.
Traceback (most recent call last):
  File "software/pylith-debug/lib/python3.12/site-packages/pylith/apps/PetscApplication.py", line 55, in onComputeNodes
    self.main(*args, **kwds)
  File "software/pylith-debug/lib/python3.12/site-packages/pylith/apps/PyLithApp.py", line 85, in main
    self.problem.initialize()
  File "software/pylith-debug/lib/python3.12/site-packages/pylith/problems/Problem.py", line 212, in initialize
    ModuleProblem.initialize(self)
  File "software/pylith-debug/lib/python3.12/site-packages/pylith/problems/problems.py", line 165, in initialize
    return _problems.Problem_initialize(self)
           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
RuntimeError: Error detected while in PETSc function.
C++ traceback (13 frames):
  [0]  software/pylith-debug/lib/libpylith.so.0(_ZN6pylith8topology10FieldQuery7queryDBEv+0x292) [0x75c8149a1782]
  [1]  software/pylith-debug/lib/libpylith.so.0(_ZN6pylith10feassemble16AuxiliaryFactory15setValuesFromDBEv+0x2ca) [0x75c8148025e4]
  [2]  software/pylith-debug/lib/libpylith.so.0(_ZN6pylith6faults6KinSrc10initializeERKNS_8topology5FieldERKNS_6scales6ScalesEPKN11spatialdata9geocoords8CoordSysE+0x51d) [0x75c8147880cd]
  [3]  software/pylith-debug/lib/libpylith.so.0(_ZN6pylith6faults16FaultCohesiveKin20createAuxiliaryFieldERKNS_8topology5FieldERKNS2_4MeshE+0xa51) [0x75c81476f305]
  [4]  software/pylith-debug/lib/libpylith.so.0(_ZN6pylith10feassemble10Integrator10initializeERKNS_8topology5FieldE+0x29f) [0x75c8147a7693]
  [5]  software/pylith-debug/lib/libpylith.so.0(_ZN6pylith10feassemble19IntegratorInterface10initializeERKNS_8topology5FieldE+0x3d6) [0x75c8147ca87c]
  [6]  software/pylith-debug/lib/libpylith.so.0(_ZN6pylith8problems7Problem10initializeEv+0x400) [0x75c81491d846]
  [7]  software/pylith-debug/lib/libpylith.so.0(_ZN6pylith8problems13TimeDependent10initializeEv+0x26f) [0x75c81492d011]
  [8]  software/pylith-debug/lib/python3.12/site-packages/pylith/problems/_problems.so(+0x15c30) [0x75c803d7ac30]
  [9]  software/pylith-debug/bin/mpinemesis(+0x15a51) [0x5717b3696a51]
  [10]  /lib/x86_64-linux-gnu/libc.so.6(+0x29d90) [0x75c827029d90]
  [11]  /lib/x86_64-linux-gnu/libc.so.6(__libc_start_main+0x80) [0x75c827029e40]
  [12]  software/pylith-debug/bin/mpinemesis(+0x5805) [0x5717b3686805]

--------------------------------------------------------------------------
MPI_ABORT was invoked on rank 0 in communicator MPI_COMM_WORLD
  Proc: [[4252,1],0]
  Errorcode: -1

NOTE: invoking MPI_ABORT causes Open MPI to kill all MPI processes.
You may or may not see output from other processes, depending on
exactly when Open MPI kills them.
--------------------------------------------------------------------------
--------------------------------------------------------------------------
prterun has exited due to process rank 0 with PID 0 on node igskci164warlng calling
"abort". This may have caused other processes in the application to be
terminated by signals sent by prterun (as reported here).
--------------------------------------------------------------------------
software/pylith-debug/bin/nemesis: mpiexec: exit 255
software/pylith-debug/bin/pylith: software/pylith-debug/bin/nemesis: exit 1

```

## Troubleshooting Strategy

With the linear interpolation we get an error about not being able to find an initiation time for a point.
This suggests there are errors in our spatial database file related to interpolation.
We examine the header and data points for errors.
We notice that our points lie along a line (data dimension is 1), but our header has `data-dim=2`.

## Resolution

```{code-block} cfg
---
caption: Correct error in `fault_slip.spatialdb`.
---
# Error
data-dim = 2
    
# Correct
data-dim = 1
```
