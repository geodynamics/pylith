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
 >> src/cig/pylith/libsrc/pylith/topology/FieldOps.cc:224:static void pylith::topology::FieldOps::checkSubfieldsExist(const pylith::string_vector &, const std::string &, const pylith::topology::Field &)
 -- error (user-input)
 -- Could not find 'lagrange_multiplier_fault' in domain solution field. Field contains: 'displacement'; the missing fields are required for interface 'FaultCohesiveKin'.
Fatal error. Calling MPI_Abort() to abort PyLith application.
Traceback (most recent call last):
  File "software/pylith-debug/lib/python3.12/site-packages/pylith/apps/PetscApplication.py", line 55, in onComputeNodes
    self.main(*args, **kwds)
  File "software/pylith-debug/lib/python3.12/site-packages/pylith/apps/PyLithApp.py", line 84, in main
    self.problem.verifyConfiguration()
  File "software/pylith-debug/lib/python3.12/site-packages/pylith/problems/Problem.py", line 206, in verifyConfiguration
    ModuleProblem.verifyConfiguration(self)
  File "software/pylith-debug/lib/python3.12/site-packages/pylith/problems/problems.py", line 162, in verifyConfiguration
    return _problems.Problem_verifyConfiguration(self)
           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
RuntimeError: Could not find 'lagrange_multiplier_fault' in domain solution field. Field contains: 'displacement'; the missing fields are required for interface 'FaultCohesiveKin'.
C++ traceback (8 frames):
  [0]  3   libpylith.0.dylib                   0x0000000108002e40 _ZN6pylith10ValueErrorCI1NS_5ErrorEERKNS_12ErrorMessageE + 36
  [1]  4   libpylith.0.dylib                   0x000000010824f548 _ZN6pylith8topology8FieldOps19checkSubfieldsExistERKNSt3__16vectorINS2_12basic_stringIcNS2_11char_traitsIcEENS2_9allocatorIcEEEENS7_IS9_EEEERKS9_RKNS0_5FieldE + 1212
  [2]  5   libpylith.0.dylib                   0x000000010803dc14 _ZNK6pylith6faults16FaultCohesiveKin19verifyConfigurationERKNS_8topology5FieldE + 648
  [3]  6   libpylith.0.dylib                   0x00000001081d8068 _ZNK6pylith8problems7Problem19verifyConfigurationEv + 1024
  [4]  7   libpylith.0.dylib                   0x00000001081e8f4c _ZNK6pylith8problems13TimeDependent19verifyConfigurationEv + 548
  [5]  8   _problems.so                        0x000000010895e218 _ZL33_wrap_Problem_verifyConfigurationP7_objectS0_ + 216
  [6]  42  mpinemesis                          0x0000000100b91d70 main + 616
  [7]  43  dyld                                0x0000000197e0eb98 start + 6076

Abort(-1) on node 0 (rank 0 in comm 0): application called MPI_Abort(MPI_COMM_WORLD, -1) - process 0
software/pylith-debug/bin/nemesis: mpiexec: exit 255
software/pylith-debug/bin/pylith: software/pylith-debug/bin/nemesis: exit 1
```

## Troubleshooting Strategy

After the Python Traceback, we see the error message on lines 64-65.
The faults check to make sure the solution field contains the necessary subfields, and the fault cannot find the `lagrange_multiplier_fault` subfield.
The easiest way to diagnose an error like this is to view the JSON file automatically generated by PyLith; it contains all of the parameters, including any defaults used.
We point our web browser to https://geodynamics.github.io/pylith_parameters/ and load the parameter file `output/step06_twofaults-parameters.json`.
In the left panel we navigate to the solution field and see that the subfields are set to `pylith.problems.SolnDisp`, so that the solution field only contains a single subfield, `displacement`.
We want the solution field to contain both `displacement` and `lagrange_multiplier_fault`.

## Resolution

```{code-block} cfg
---
caption: Correct error in `step06_twofaults.cfg`.
---
[pylithapp.problem]
solution = pylith.problems.SolnDispLagrange
