# Step 6: Error 6

## Error Message

```{code-block} console
---
caption: Error message 6 when running Step 6.
linenos: True
emphasize-lines: 55-59, 72-75
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
 >> src/cig/pylith/libsrc/pylith/topology/FieldQuery.cc:473:static void pylith::topology::_FieldQuery::findQueryIndices(pylith::topology::FieldQuery::DBQueryContext*, const pylith::string_vector&)
 -- error (user-input)
 -- Could not find value 'final_slip_opening' in spatial database 'Fault rupture for main fault'. Available values are:
  final-slip-left-lateral
  final-slip-opening
  initiation-time

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
RuntimeError: Could not find value 'final_slip_opening' in spatial database 'Fault rupture for main fault'. Available values are:
  final-slip-left-lateral
  final-slip-opening
  initiation-time

C++ traceback (14 frames):
  [0]  software/pylith-debug/lib/libpylith.so.0(_ZN6pylith8topology11_FieldQuery16findQueryIndicesEPNS0_10FieldQuery14DBQueryContextERKSt6vectorINSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEESaISB_EE+0x60a) [0x7288b51a486e]
  [1]  software/pylith-debug/lib/libpylith.so.0(_ZN6pylith8topology10FieldQuery6openDBEPN11spatialdata9spatialdb9SpatialDBEd+0x650) [0x7288b51a0f44]
  [2]  software/pylith-debug/lib/libpylith.so.0(_ZN6pylith10feassemble16AuxiliaryFactory15setValuesFromDBEv+0x2b4) [0x7288b50025ce]
  [3]  software/pylith-debug/lib/libpylith.so.0(_ZN6pylith6faults6KinSrc10initializeERKNS_8topology5FieldERKNS_6scales6ScalesEPKN11spatialdata9geocoords8CoordSysE+0x51d) [0x7288b4f880cd]
  [4]  software/pylith-debug/lib/libpylith.so.0(_ZN6pylith6faults16FaultCohesiveKin20createAuxiliaryFieldERKNS_8topology5FieldERKNS2_4MeshE+0xa51) [0x7288b4f6f305]
  [5]  software/pylith-debug/lib/libpylith.so.0(_ZN6pylith10feassemble10Integrator10initializeERKNS_8topology5FieldE+0x29f) [0x7288b4fa7693]
  [6]  software/pylith-debug/lib/libpylith.so.0(_ZN6pylith10feassemble19IntegratorInterface10initializeERKNS_8topology5FieldE+0x3d6) [0x7288b4fca87c]
  [7]  software/pylith-debug/lib/libpylith.so.0(_ZN6pylith8problems7Problem10initializeEv+0x400) [0x7288b511d846]
  [8]  software/pylith-debug/lib/libpylith.so.0(_ZN6pylith8problems13TimeDependent10initializeEv+0x26f) [0x7288b512d011]
  [9]  software/pylith-debug/lib/python3.12/site-packages/pylith/problems/_problems.so(+0x15c30) [0x7288ae542c30]
  [10]  software/pylith-debug/bin/mpinemesis(+0x15a51) [0x5c2fcc5aaa51]
  [11]  /lib/x86_64-linux-gnu/libc.so.6(+0x29d90) [0x7288d1829d90]
  [12]  /lib/x86_64-linux-gnu/libc.so.6(__libc_start_main+0x80) [0x7288d1829e40]
  [13]  software/pylith-debug/bin/mpinemesis(+0x5805) [0x5c2fcc59a805]

--------------------------------------------------------------------------
MPI_ABORT was invoked on rank 0 in communicator MPI_COMM_WORLD
  Proc: [[30692,1],0]
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
