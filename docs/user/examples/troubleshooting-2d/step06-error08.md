# Step 6: Error 8

## Error Message

```{code-block} console
---
caption: Error message 8 when running Step 6.
linenos: True
emphasize-lines: 58-59
---
$ pylith step06_twofaults.cfg

 >> software/pylith-debug/lib/python3.12/site-packages/pylith/apps/PyLithApp.py:76:main
 -- info (application-flow)
 -- Running on 1 process(es).
 >> src/cig/pylith/libsrc/pylith/utils/PetscOptions.cc:251:static void pylith::utils::_PetscOptions::write(pythia::journal::info_t &, const char *, const PetscOptions &)
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

 >> src/cig/pylith/libsrc/pylith/meshio/MeshIOPetsc.cc:205:virtual void pylith::meshio::MeshIOPetsc::_read()
 -- info (application-flow)
 -- Component 'meshiopetsc.reader': Reading finite-element mesh from 'mesh_tri.msh'.
 >> src/cig/pylith/libsrc/pylith/meshio/MeshIO.cc:76:void pylith::meshio::MeshIO::read(pylith::topology::Mesh *, const bool)
 -- info (application-flow)
 -- Component 'meshiopetsc.reader': Domain bounding box:
    (-100000, 100000)
    (-100000, 0)
 >> src/cig/pylith/libsrc/pylith/initializers/MeshInsertInterfaces.cc:51:virtual pylith::topology::Mesh *pylith::initializers::MeshInsertInterfaces::run(pylith::topology::Mesh *, const pylith::problems::Problem &)
 -- info (application-flow)
 -- Inserting cohesive cells.
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
[0]PETSC ERROR: PETSc Development Git Revision: v3.25.1-142-g990e4f00326 Git Date: 2026-05-15 01:03:36 -0400
[0]PETSC ERROR: software/pylith-debug/bin/mpinemesis with 1 MPI process(es) and PETSC_ARCH arch-pylith-debug on IGSKCI164LM006 by baagaard Thu Jun  4 12:15:03 2026
[0]PETSC ERROR: Configure options: --PETSC_ARCH=arch-pylith-debug --with-debugging=1 --with-clanguage=c --with-mpi-compilers=1 --with-shared-libraries=1 --with-64-bit-points=1 --with-large-file-io=1 --with-lgrind=0 --download-parmetis=1 --download-metis=1 --download-triangle --download-ml=1 --download-superlu=1 --with-fc=0 --download-f2cblaslapack --with-hdf5=1 --with-hdf5-include=software/pylith-debug/include --with-hdf5-lib=software/pylith-debug/lib/libhdf5.dylib --with-zlib=1 CFLAGS+=-g
[0]PETSC ERROR: #1 static PetscErrorCode pylith::topology::FieldQuery::queryDBPointFn(PylithInt, PylithReal, const PylithReal *, PylithInt, PylithScalar *, void *)() at src/cig/pylith/libsrc/pylith/topology/FieldQuery.cc:320
[0]PETSC ERROR: #2 DMProjectPoint_Func_Private() at software/unix/petsc-dev/src/dm/impls/plex/plexproject.c:143
[0]PETSC ERROR: #3 DMProjectPoint_Private() at software/unix/petsc-dev/src/dm/impls/plex/plexproject.c:545
[0]PETSC ERROR: #4 DMProjectLocal_Generic_Plex() at software/unix/petsc-dev/src/dm/impls/plex/plexproject.c:1083
[0]PETSC ERROR: #5 DMProjectFunctionLocal_Plex() at software/unix/petsc-dev/src/dm/impls/plex/plexproject.c:1114
[0]PETSC ERROR: #6 DMProjectFunctionLocal() at software/unix/petsc-dev/src/dm/interface/dm.c:8361
[0]PETSC ERROR: #7 void pylith::topology::FieldQuery::queryDB()() at src/cig/pylith/libsrc/pylith/topology/FieldQuery.cc:228
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
  [0]  3   libpylith.0.dylib                   0x000000010653c9fc _ZN6pylith13ExternalErrorCI1NS_5ErrorEERKNS_12ErrorMessageE + 36
  [1]  4   libpylith.0.dylib                   0x0000000106792c28 _ZN6pylith8topology10FieldQuery7queryDBEv + 592
  [2]  5   libpylith.0.dylib                   0x000000010660309c _ZN6pylith10feassemble16AuxiliaryFactory15setValuesFromDBEv + 608
  [3]  6   libpylith.0.dylib                   0x000000010658ed48 _ZN6pylith6faults6KinSrc10initializeERKNS_8topology5FieldERKNS_6scales6ScalesEPKN11spatialdata9geocoords8CoordSysE + 1180
  [4]  7   libpylith.0.dylib                   0x0000000106576f64 _ZN6pylith6faults16FaultCohesiveKin20createAuxiliaryFieldERKNS_8topology5FieldERKNS2_4MeshE + 2764
  [5]  8   libpylith.0.dylib                   0x00000001065ac9bc _ZN6pylith10feassemble10Integrator10initializeERKNS_8topology5FieldE + 560
  [6]  9   libpylith.0.dylib                   0x00000001065cd944 _ZN6pylith10feassemble19IntegratorInterface10initializeERKNS_8topology5FieldE + 924
  [7]  10  libpylith.0.dylib                   0x0000000106711014 _ZN6pylith8problems7Problem10initializeEv + 936
  [8]  11  libpylith.0.dylib                   0x0000000106721610 _ZN6pylith8problems13TimeDependent10initializeEv + 652
  [9]  12  _problems.so                        0x00000001096e23f0 _ZL24_wrap_Problem_initializeP7_objectS0_ + 216
  [10]  46  mpinemesis                          0x0000000102addd70 main + 616
  [11]  47  dyld                                0x0000000197e0eb98 start + 6076

Abort(-1) on node 0 (rank 0 in comm 0): application called MPI_Abort(MPI_COMM_WORLD, -1) - process 0
software/pylith-debug/bin/nemesis: mpiexec: exit 255
software/pylith-debug/bin/pylith: software/pylith-debug/bin/nemesis: exit 1
```

## Troubleshooting Strategy

We still get an error about not being able to find an initiation time for a point.
This suggests there are still one or more errors in our spatial database file related to interpolation.
We examine the header and data points for errors.
We notice that our deepest point has a y coordinate of -25 km, but PyLith is looking for values at a point with a y coordinate of -27.621 km.
We need to add an additional point to our spatial database.
This explains why we had `num-locs=4` when we started!

## Resolution

```{code-block} cfg
---
caption: Correct error in `fault_slip.spatialdb`.
---
num-locs = 4
...
0.0   99.0     -2.0       0.0   0.0
0.0  -20.0     -2.0       0.0   0.0
0.0  -25.0      0.0       0.0   0.0
0.0  -99.0      0.0       0.0   0.0
```
