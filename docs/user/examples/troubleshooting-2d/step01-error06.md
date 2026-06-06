# Step 1: Successful run

```{code-block} pyrejournal
---
caption: Output when running Step 1.
linenos: True
---
$ pylith step01a_gravity.cfg

# Output
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
 >> src/cig/pylith/libsrc/pylith/problems/TimeDependent.cc:473:void pylith::problems::TimeDependent::solve()
 -- info (application-flow)
 -- Component 'timedependent.problem': Solving equations.
0 TS dt 0.001 time 0.
    0 SNES Function norm 3.007829881319e+03
      Linear solve converged due to CONVERGED_ATOL iterations 11
    1 SNES Function norm 3.524261285349e-07
    Nonlinear solve converged due to CONVERGED_FNORM_ABS iterations 1
1 TS dt 0.001 time 0.001
 >> software/pylith-debug/lib/python3.12/site-packages/pylith/problems/Problem.py:222:finalize
 -- info (application-flow)
 -- Finalizing problem.
```

Our simulation now runs without errors and the output looks correct.
