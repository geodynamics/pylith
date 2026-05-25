# Step 6: Error 3

## Error Message

```{code-block} console
---
caption: Error message 3 when running Step 6.
linenos: True
emphasize-lines: 41-42, 73-74
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
[0]PETSC ERROR: --------------------- Error Message --------------------------------------------------------------
[0]PETSC ERROR: Invalid argument
[0]PETSC ERROR: Edge has 2 coordinates != 4
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
[0]PETSC ERROR:   Option left: name:-viewer_hdf5_collective (no value) source: code
[0]PETSC ERROR: See https://petsc.org/release/faq/ for trouble shooting.
[0]PETSC ERROR: PETSc Development Git Revision: v3.25.1-168-g70613f6caab Git Date: 2026-05-21 15:15:28 +0000
[0]PETSC ERROR: software/pylith-debug/bin/mpinemesis with 1 MPI process(es) and PETSC_ARCH arch-pylith-debug on igskci164warlng.gs.doi.net by baagaard Tue May 26 11:56:50 2026
[0]PETSC ERROR: Configure options: --PETSC_ARCH=arch-pylith-debug --with-debugging=1 --with-clanguage=c --with-mpi-compilers=1 --with-shared-libraries=1 --with-64-bit-points=1 --with-large-file-io=1 --with-lgrind=0 --download-parmetis=1 --download-metis=1 --download-triangle --download-ml=1 --download-superlu=1 --with-fc=0 --download-f2cblaslapack --with-hdf5=1 --with-hdf5-dir=software/pylith-debug --with-zlib=1 CFLAGS+=-g
[0]PETSC ERROR: #1 DMPlexComputeGeometryFVM_1D_Internal() at /software/baagaard/petsc-dev/src/dm/impls/plex/plexgeometry.c:2692
[0]PETSC ERROR: #2 DMPlexComputeCellGeometryFVM() at /software/baagaard/petsc-dev/src/dm/impls/plex/plexgeometry.c:2935
[0]PETSC ERROR: #3 DMPlexLabelCohesiveCheck() at /software/baagaard/petsc-dev/src/dm/impls/plex/plexsubmesh.c:2992
[0]PETSC ERROR: #4 static void pylith::faults::TopologyOps::updateCohesiveLabel(const pylith::topology::Mesh*, const char*, int)() at src/cig/pylith/libsrc/pylith/faults/TopologyOps.cc:86
Fatal error. Calling MPI_Abort() to abort PyLith application.
Traceback (most recent call last):
  File "software/pylith-debug/lib/python3.12/site-packages/pylith/apps/PetscApplication.py", line 55, in onComputeNodes
    self.main(*args, **kwds)
  File "software/pylith-debug/lib/python3.12/site-packages/pylith/apps/PyLithApp.py", line 83, in main
    self.problem.preinitialize()
  File "software/pylith-debug/lib/python3.12/site-packages/pylith/problems/TimeDependent.py", line 120, in preinitialize
    Problem.preinitialize(self)
  File "software/pylith-debug/lib/python3.12/site-packages/pylith/problems/Problem.py", line 190, in preinitialize
    self.mesh = self.meshInitializer.runPhases(self)
                ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "software/pylith-debug/lib/python3.12/site-packages/pylith/initializers/initializers.py", line 76, in runPhases
    return _initializers.Initializer_runPhases(self, problem)
           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
RuntimeError: Error detected while in PETSc function.
Error occurred while transforming topology to create cohesive cells for fault 'fault'.

C++ traceback (9 frames):
  [0]  software/pylith-debug/lib/libpylith.so.0(_ZN6pylith6faults11TopologyOps19updateCohesiveLabelEPKNS_8topology4MeshEPKci+0x22fb) [0x7b82a379c16f]
  [1]  software/pylith-debug/lib/libpylith.so.0(_ZN6pylith6faults13FaultCohesive17transformTopologyEPNS_8topology4MeshE+0x854) [0x7b82a376336a]
  [2]  software/pylith-debug/lib/libpylith.so.0(_ZN6pylith12initializers20MeshInsertInterfaces3runEPNS_8topology4MeshERKNS_8problems7ProblemE+0xa03) [0x7b82a380e58b]
  [3]  software/pylith-debug/lib/libpylith.so.0(_ZN6pylith12initializers11Initializer9runPhasesERKNS_8problems7ProblemE+0x2e2) [0x7b82a3806736]
  [4]  software/pylith-debug/lib/python3.12/site-packages/pylith/initializers/_initializers.so(+0xb409) [0x7b82897cf409]
  [5]  software/pylith-debug/bin/mpinemesis(+0x15a51) [0x5a53947afa51]
  [6]  /lib/x86_64-linux-gnu/libc.so.6(+0x29d90) [0x7b82c0029d90]
  [7]  /lib/x86_64-linux-gnu/libc.so.6(__libc_start_main+0x80) [0x7b82c0029e40]
  [8]  software/pylith-debug/bin/mpinemesis(+0x5805) [0x5a539479f805]

--------------------------------------------------------------------------
MPI_ABORT was invoked on rank 0 in communicator MPI_COMM_WORLD
  Proc: [[26341,1],0]
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

We no longer have errors during the problem configuration.
Now we have errors while creating the cohesive cells.
PETSc reports an error on lines 64-65.
The message indicates an error with the resulting topology after creating cohesive cells for fault 'fault'.
In `step06_twofaults.cfg` we see that the cohesive cells are created for fault 'splay' and then for fault 'fault'.
Fault 'fault' is the main fault and should be created first, followed by fault 'splay'.

## Resolution

```{code-block} cfg
---
caption: Correct error in `step06_twofaults.cfg`.
---
# Errpr
[pylithapp.problem]
interfaces = [splay, fault]

# Correct
[pylithapp.problem]
interfaces = [fault, splay]
```
