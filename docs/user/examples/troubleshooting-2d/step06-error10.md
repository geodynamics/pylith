# Step 6: Successful run

```{code-block} console
---
caption: Output when running Step 6.
linenos: True
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
 >> src/cig/pylith/libsrc/pylith/problems/TimeDependent.cc:473:void pylith::problems::TimeDependent::solve()
 -- info (application-flow)
 -- Component 'timedependent.problem': Solving equations.
0 TS dt 0.2 time -0.2
    0 SNES Function norm 7.713717485060e+00
      Linear solve converged due to CONVERGED_ATOL iterations 14
    1 SNES Function norm 5.241642052618e-08
    Nonlinear solve converged due to CONVERGED_FNORM_ABS iterations 1
1 TS dt 0.2 time 0.
    0 SNES Function norm 5.241642052618e-08
    Nonlinear solve converged due to CONVERGED_FNORM_ABS iterations 0
2 TS dt 0.2 time 0.2
    0 SNES Function norm 3.245094761956e+00
      Linear solve converged due to CONVERGED_ATOL iterations 13
    1 SNES Function norm 7.651348491638e-08
    Nonlinear solve converged due to CONVERGED_FNORM_ABS iterations 1
3 TS dt 0.2 time 0.4
 >> software/pylith-debug/lib/python3.12/site-packages/pylith/problems/Problem.py:222:finalize
 -- info (application-flow)
 -- Finalizing problem.
```

Our simulation now runs without errors and the output looks correct.
