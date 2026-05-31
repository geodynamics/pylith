# Step 6: Error 4

## Error Message

```{code-block} console
---
caption: Error message 4 when running Step 6.
linenos: True
emphasize-lines: 44-45, 57
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
 >> src/cig/pylith/libsrc/pylith/topology/FieldOps.cc:223:static void pylith::topology::FieldOps::checkSubfieldsExist(const pylith::string_vector&, const std::string&, const pylith::topology::Field&)
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
C++ traceback (9 frames):
  [0]  software/pylith-debug/lib/libpylith.so.0(_ZN6pylith8topology8FieldOps19checkSubfieldsExistERKSt6vectorINSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEESaIS8_EERKS8_RKNS0_5FieldE+0x480) [0x7eb3a43984c0]
  [1]  software/pylith-debug/lib/libpylith.so.0(_ZNK6pylith6faults16FaultCohesiveKin19verifyConfigurationERKNS_8topology5FieldE+0x37b) [0x7eb3a416dfff]
  [2]  software/pylith-debug/lib/libpylith.so.0(_ZNK6pylith8problems7Problem19verifyConfigurationEv+0x474) [0x7eb3a431cf68]
  [3]  software/pylith-debug/lib/libpylith.so.0(_ZNK6pylith8problems13TimeDependent19verifyConfigurationEv+0x253) [0x7eb3a432c8e5]
  [4]  software/pylith-debug/lib/python3.12/site-packages/pylith/problems/_problems.so(+0x15af4) [0x7eb39de12af4]
  [5]  software/pylith-debug/bin/mpinemesis(+0x15a51) [0x5bd4a7faaa51]
  [6]  /lib/x86_64-linux-gnu/libc.so.6(+0x29d90) [0x7eb3c0a29d90]
  [7]  /lib/x86_64-linux-gnu/libc.so.6(__libc_start_main+0x80) [0x7eb3c0a29e40]
  [8]  software/pylith-debug/bin/mpinemesis(+0x5805) [0x5bd4a7f9a805]

--------------------------------------------------------------------------
MPI_ABORT was invoked on rank 0 in communicator MPI_COMM_WORLD
  Proc: [[52992,1],0]
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
