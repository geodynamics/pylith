# Step 6: Error 5

## Error Message

```{code-block} console
---
caption: Error message 5 when running Step 6.
linenos: True
emphasize-lines: 58-62,75-78
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
 >> src/cig/pylith/libsrc/pylith/topology/FieldQuery.cc:473:static void pylith::topology::_FieldQuery::findQueryIndices(FieldQuery::DBQueryContext *, const pylith::string_vector &)
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

C++ traceback (13 frames):
  [0]  3   libpylith.0.dylib                   0x00000001041c2e40 _ZN6pylith10ValueErrorCI1NS_5ErrorEERKNS_12ErrorMessageE + 36
  [1]  4   libpylith.0.dylib                   0x000000010441a864 _ZN6pylith8topology11_FieldQuery16findQueryIndicesEPNS0_10FieldQuery14DBQueryContextERKNSt3__16vectorINS5_12basic_stringIcNS5_11char_traitsIcEENS5_9allocatorIcEEEENSA_ISC_EEEE + 1624
  [2]  5   libpylith.0.dylib                   0x0000000104418ae4 _ZN6pylith8topology10FieldQuery6openDBEPN11spatialdata9spatialdb9SpatialDBEd + 1720
  [3]  6   libpylith.0.dylib                   0x000000010428b090 _ZN6pylith10feassemble16AuxiliaryFactory15setValuesFromDBEv + 596
  [4]  7   libpylith.0.dylib                   0x0000000104216d48 _ZN6pylith6faults6KinSrc10initializeERKNS_8topology5FieldERKNS_6scales6ScalesEPKN11spatialdata9geocoords8CoordSysE + 1180
  [5]  8   libpylith.0.dylib                   0x00000001041fef64 _ZN6pylith6faults16FaultCohesiveKin20createAuxiliaryFieldERKNS_8topology5FieldERKNS2_4MeshE + 2764
  [6]  9   libpylith.0.dylib                   0x00000001042349bc _ZN6pylith10feassemble10Integrator10initializeERKNS_8topology5FieldE + 560
  [7]  10  libpylith.0.dylib                   0x0000000104255944 _ZN6pylith10feassemble19IntegratorInterface10initializeERKNS_8topology5FieldE + 924
  [8]  11  libpylith.0.dylib                   0x0000000104399014 _ZN6pylith8problems7Problem10initializeEv + 936
  [9]  12  libpylith.0.dylib                   0x00000001043a9610 _ZN6pylith8problems13TimeDependent10initializeEv + 652
  [10]  13  _problems.so                        0x00000001079ce3f0 _ZL24_wrap_Problem_initializeP7_objectS0_ + 216
  [11]  47  mpinemesis                          0x00000001006e9d70 main + 616
  [12]  48  dyld                                0x0000000197e0eb98 start + 6076

Abort(-1) on node 0 (rank 0 in comm 0): application called MPI_Abort(MPI_COMM_WORLD, -1) - process 0
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
