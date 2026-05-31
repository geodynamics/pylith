# Step 1: Error 5

## Error Message

```{code-block} console
---
caption: Error message 5 when running Step 1.
linenos: True
emphasize-lines: 52-53
---
$ pylith step01a_gravity.cfg

 >> software/pylith-debug/lib/python3.12/site-packages/pylith/apps/PyLithApp.py:76:main
 -- info (application-flow)
 -- Running on 1 process(es).
 >> src/cig/pylith/libsrc/pylith/utils/PetscOptions.cc:251:static void pylith::utils::_PetscOptions::write(pythia::journal::info_t&, const char*, const pylith::utils::PetscOptions&)
 -- info (application-flow)
 -- Setting PETSc options:
    ksp_atol = 1.0e-7
    ksp_converged_reason = true
    ksp_error_if_not_converged = true
    ksp_gmres_restart = 100
    ksp_guess_pod_size = 8
    ksp_guess_type = pod
    ksp_rtol = 1.0e-14
    mg_fine_ksp_max_it = 5
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
    Length scale: 5000*m
    Displacement scale: 1*m
    Time scale: 3.15576e+09*s
    Rigidity scale: 1e+10*m**-1*kg*s**-2
    Temperature scale: 1*K
 >> src/cig/pylith/libsrc/pylith/problems/TimeDependent.cc:342:virtual void pylith::problems::TimeDependent::initialize()
 -- info (application-flow)
 -- Component 'timedependent.problem': Initializing problem.
[0]PETSC ERROR: --------------------- Error Message --------------------------------------------------------------
[0]PETSC ERROR: Error in external library
[0]PETSC ERROR: Error converting spatial database values for gravitational_acceleration at (  -92285.3  -72772.7) in spatial database 'Gravity field'. Found near zero magnitude (0) for gravity field vector (  0  0).
[0]PETSC ERROR: WARNING! There are unused option(s) set! Could be the program crashed before usage or a spelling mistake, etc!
[0]PETSC ERROR:   Option left: name:-ksp_atol value: 1.0e-7 source: code
[0]PETSC ERROR:   Option left: name:-ksp_converged_reason (no value) source: code
[0]PETSC ERROR:   Option left: name:-ksp_error_if_not_converged (no value) source: code
[0]PETSC ERROR:   Option left: name:-ksp_gmres_restart value: 100 source: code
[0]PETSC ERROR:   Option left: name:-ksp_guess_pod_size value: 8 source: code
[0]PETSC ERROR:   Option left: name:-ksp_guess_type value: pod source: code
[0]PETSC ERROR:   Option left: name:-ksp_rtol value: 1.0e-14 source: code
[0]PETSC ERROR:   Option left: name:-mg_fine_ksp_max_it value: 5 source: code
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
[0]PETSC ERROR:   Option left: name:-viewer_hdf5_collective (no value) source: code
[0]PETSC ERROR: See https://petsc.org/release/faq/ for trouble shooting.
[0]PETSC ERROR: PETSc Development Git Revision: v3.25.1-168-g70613f6caab Git Date: 2026-05-21 15:15:28 +0000
[0]PETSC ERROR: software/pylith-debug/bin/mpinemesis with 1 MPI process(es) and PETSC_ARCH arch-pylith-debug on igskci164warlng.gs.doi.net by baagaard Tue May 26 11:31:47 2026
[0]PETSC ERROR: Configure options: --PETSC_ARCH=arch-pylith-debug --with-debugging=1 --with-clanguage=c --with-mpi-compilers=1 --with-shared-libraries=1 --with-64-bit-points=1 --with-large-file-io=1 --with-lgrind=0 --download-parmetis=1 --download-metis=1 --download-triangle --download-ml=1 --download-superlu=1 --with-fc=0 --download-f2cblaslapack --with-hdf5=1 --with-hdf5-dir=software/pylith-debug --with-zlib=1 CFLAGS+=-g
[0]PETSC ERROR: #1 static PetscErrorCode pylith::topology::FieldQuery::queryDBPointFn(PylithInt, PylithReal, const PylithReal*, PylithInt, PylithScalar*, void*)() at src/cig/pylith/libsrc/pylith/topology/FieldQuery.cc:334
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
C++ traceback (12 frames):
  [0]  software/pylith-debug/lib/libpylith.so.0(_ZN6pylith8topology10FieldQuery7queryDBEv+0x292) [0x7a7abcda1782]
  [1]  software/pylith-debug/lib/libpylith.so.0(_ZN6pylith10feassemble16AuxiliaryFactory15setValuesFromDBEv+0x2ca) [0x7a7abcc025e4]
  [2]  software/pylith-debug/lib/libpylith.so.0(_ZN6pylith9materials10Elasticity20createAuxiliaryFieldERKNS_8topology5FieldERKNS2_4MeshE+0x561) [0x7a7abcc154e1]
  [3]  software/pylith-debug/lib/libpylith.so.0(_ZN6pylith10feassemble10Integrator10initializeERKNS_8topology5FieldE+0x29f) [0x7a7abcba7693]
  [4]  software/pylith-debug/lib/libpylith.so.0(_ZN6pylith10feassemble16IntegratorDomain10initializeERKNS_8topology5FieldE+0x3b1) [0x7a7abcbb1def]
  [5]  software/pylith-debug/lib/libpylith.so.0(_ZN6pylith8problems7Problem10initializeEv+0x400) [0x7a7abcd1d846]
  [6]  software/pylith-debug/lib/libpylith.so.0(_ZN6pylith8problems13TimeDependent10initializeEv+0x26f) [0x7a7abcd2d011]
  [7]  software/pylith-debug/lib/python3.12/site-packages/pylith/problems/_problems.so(+0x15c30) [0x7a7ab67c4c30]
  [8]  software/pylith-debug/bin/mpinemesis(+0x15a51) [0x5fc5abcf9a51]
  [9]  /lib/x86_64-linux-gnu/libc.so.6(+0x29d90) [0x7a7ad9429d90]
  [10]  /lib/x86_64-linux-gnu/libc.so.6(__libc_start_main+0x80) [0x7a7ad9429e40]
  [11]  software/pylith-debug/bin/mpinemesis(+0x5805) [0x5fc5abce9805]

--------------------------------------------------------------------------
MPI_ABORT was invoked on rank 0 in communicator MPI_COMM_WORLD
  Proc: [[30480,1],0]
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
