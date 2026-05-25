# Step 1: Error 4

## Error Message

```{code-block} console
---
caption: Error message 4 when running Step 1.
linenos: True
emphasize-lines: 41-42, 54
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
 >> src/cig/pylith/libsrc/pylith/topology/MeshOps.cc:802:static void pylith::topology::MeshOps::checkMaterialLabels(const pylith::topology::Mesh&, pylith::int_array&)
 -- error (user-input)
 -- Material label_value '3' for cell '609' does not match the label_value of any materials or interfaces.
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
RuntimeError: Material label_value '3' for cell '609' does not match the label_value of any materials or interfaces.
C++ traceback (9 frames):
  [0]  software/pylith-debug/lib/libpylith.so.0(_ZN6pylith8topology7MeshOps19checkMaterialLabelsERKNS0_4MeshERSt8valarrayIiE+0x7b9) [0x74544c5780c7]
  [1]  software/pylith-debug/lib/libpylith.so.0(_ZNK6pylith8problems7Problem20_checkMaterialLabelsEv+0x4f6) [0x74544c51e69e]
  [2]  software/pylith-debug/lib/libpylith.so.0(_ZNK6pylith8problems7Problem19verifyConfigurationEv+0x552) [0x74544c51d046]
  [3]  software/pylith-debug/lib/libpylith.so.0(_ZNK6pylith8problems13TimeDependent19verifyConfigurationEv+0x253) [0x74544c52c8e5]
  [4]  software/pylith-debug/lib/python3.12/site-packages/pylith/problems/_problems.so(+0x15af4) [0x745445fd8af4]
  [5]  software/pylith-debug/bin/mpinemesis(+0x15a51) [0x612c79077a51]
  [6]  /lib/x86_64-linux-gnu/libc.so.6(+0x29d90) [0x745468c29d90]
  [7]  /lib/x86_64-linux-gnu/libc.so.6(__libc_start_main+0x80) [0x745468c29e40]
  [8]  software/pylith-debug/bin/mpinemesis(+0x5805) [0x612c79067805]

--------------------------------------------------------------------------
MPI_ABORT was invoked on rank 0 in communicator MPI_COMM_WORLD
  Proc: [[1354,1],0]
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

The output shows a Python Traceback and then the error message on Line 52.
The error indicates the finite-element mesh file contains a cell with a label value of 3, but the parameter files do not have a material with a label value of 3.
We examine the `pylithapp.problem.materials` sections of `pylithapp.cfg` and see that the label values are 0, 1, and 2 in the parameter file rather than 1, 2, and 3 that we set in the Gmsh Python script.

## Resolution

```{code-block} cfg
---
caption: Correct error in `pylithapp.cfg`.
---
[pylithapp.problem.materials.slab]
label_value = 1
...

[pylithapp.problem.materials.crust]
label_value = 2
...

[pylithapp.problem.materials.wedge]
label_value = 3
...
```
