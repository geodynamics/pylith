# Step 1: Gravitational Body Forces

This example involves a static simulation that solves for the deformation from loading by gravitational body forces.
{numref}`fig:example:reverse:2d:step01:diagram` shows the boundary conditions on the domain.

:::{figure-md} fig:example:reverse:2d:step01:diagram
<img src="figs/step01-diagram.*" alt="" scale="75%">

We apply roller boundary conditions on the lateral sides and bottom of the domain.
:::

% Metadata extracted from parameter files.
```{include} step01_gravity-synopsis.md
```

## Simulation parameters

The parameters specific to this example are in `step01_gravity.cfg`.
These include:

* `pylithapp.metadata` Metadata for this simulation. Even when the author and version are the same for all simulations in a directory, we prefer to keep that metadata in each simulation file as a reminder to keep it up-to-date for each simulation.
* `pylithapp` Parameters defining where to write the output.
* `pylithapp.problem` Parameters for specifying the gravitational body forces and adjusting the basis order.

```{code-block} console
---
caption: Run Step 1 simulation
---
$ pylith step01_gravity.cfg

# The output should look something like the following.
 >> /Users/baagaard/software/unix/py39-venv/pylith-debug/lib/python3.9/site-packages/pylith/meshio/MeshIOObj.py:44:read
 -- meshiopetsc(info)
 -- Reading finite-element mesh
 >> /Users/baagaard/src/cig/pylith/libsrc/pylith/meshio/MeshIO.cc:94:void pylith::meshio::MeshIO::read(topology::Mesh *)
 -- meshiopetsc(info)
 -- Component 'reader': Domain bounding box:
    (-100000, 100000)
    (-100000, 0)
 >> /Users/baagaard/software/unix/py39-venv/pylith-debug/lib/python3.9/site-packages/pylith/problems/Problem.py:116:preinitialize
 -- timedependent(info)
 -- Performing minimal initialization before verifying configuration.
 >> /Users/baagaard/software/unix/py39-venv/pylith-debug/lib/python3.9/site-packages/pylith/problems/Solution.py:44:preinitialize
 -- solution(info)
 -- Performing minimal initialization of solution.
 >> /Users/baagaard/software/unix/py39-venv/pylith-debug/lib/python3.9/site-packages/pylith/materials/RheologyElasticity.py:41:preinitialize
 -- isotropiclinearelasticity(info)
 -- Performing minimal initialization of elasticity rheology 'bulk_rheology'.
 >> /Users/baagaard/software/unix/py39-venv/pylith-debug/lib/python3.9/site-packages/pylith/materials/RheologyElasticity.py:41:preinitialize
 -- isotropiclinearelasticity(info)
 -- Performing minimal initialization of elasticity rheology 'bulk_rheology'.
 >> /Users/baagaard/software/unix/py39-venv/pylith-debug/lib/python3.9/site-packages/pylith/materials/RheologyElasticity.py:41:preinitialize
 -- isotropiclinearelasticity(info)
 -- Performing minimal initialization of elasticity rheology 'bulk_rheology'.
 >> /Users/baagaard/software/unix/py39-venv/pylith-debug/lib/python3.9/site-packages/pylith/bc/DirichletTimeDependent.py:92:preinitialize
 -- dirichlettimedependent(info)
 -- Performing minimal initialization of time-dependent Dirichlet boundary condition 'bc_xneg'.
 >> /Users/baagaard/software/unix/py39-venv/pylith-debug/lib/python3.9/site-packages/pylith/bc/DirichletTimeDependent.py:92:preinitialize
 -- dirichlettimedependent(info)
 -- Performing minimal initialization of time-dependent Dirichlet boundary condition 'bc_xpos'.
 >> /Users/baagaard/software/unix/py39-venv/pylith-debug/lib/python3.9/site-packages/pylith/bc/DirichletTimeDependent.py:92:preinitialize
 -- dirichlettimedependent(info)
 -- Performing minimal initialization of time-dependent Dirichlet boundary condition 'bc_yneg'.
 >> /Users/baagaard/software/unix/py39-venv/pylith-debug/lib/python3.9/site-packages/pylith/problems/Problem.py:175:verifyConfiguration
 -- timedependent(info)
 -- Verifying compatibility of problem configuration.
 >> /Users/baagaard/software/unix/py39-venv/pylith-debug/lib/python3.9/site-packages/pylith/problems/Problem.py:221:_printInfo
 -- timedependent(info)
 -- Scales for nondimensionalization:
    Length scale: 1000*m
    Time scale: 3.15576e+09*s
    Pressure scale: 3e+10*m**-1*kg*s**-2
    Density scale: 2.98765e+23*m**-3*kg
    Temperature scale: 1*K
 >> /Users/baagaard/software/unix/py39-venv/pylith-debug/lib/python3.9/site-packages/pylith/problems/Problem.py:186:initialize
 -- timedependent(info)
 -- Initializing timedependent problem with quasistatic formulation.
 >> /Users/baagaard/src/cig/pylith/libsrc/pylith/utils/PetscOptions.cc:235:static void pylith::utils::_PetscOptions::write(pythia::journal::info_t &, const char *, const pylith::utils::PetscOptions &)
 -- petscoptions(info)
 -- Setting PETSc options:
ksp_atol = 1.0e-12
ksp_converged_reason = true
ksp_error_if_not_converged = true
ksp_rtol = 1.0e-12
pc_type = lu
snes_atol = 1.0e-9
snes_converged_reason = true
snes_error_if_not_converged = true
snes_monitor = true
snes_rtol = 1.0e-12
ts_error_if_step_fails = true
ts_monitor = true
ts_type = beuler

 >> /Users/baagaard/software/unix/py39-venv/pylith-debug/lib/python3.9/site-packages/pylith/problems/TimeDependent.py:139:run
 -- timedependent(info)
 -- Solving problem.
0 TS dt 0.01 time 0.
    0 SNES Function norm 2.873918352757e-01 
    Linear solve converged due to CONVERGED_RTOL iterations 1
    1 SNES Function norm 3.025686251687e-13 
  Nonlinear solve converged due to CONVERGED_FNORM_ABS iterations 1
1 TS dt 0.01 time 0.01
 >> /Users/baagaard/software/unix/py39-venv/pylith-debug/lib/python3.9/site-packages/pylith/problems/Problem.py:201:finalize
 -- timedependent(info)
 -- Finalizing problem.
```

At the beginning of the output written to the terminal, we see that PyLith is reading the mesh using the `MeshIOPetsc` reader and that it found the domain to extend from -100,000 m to +100,000 m in the x direction and from -100,000 m to 0 in the y direction.
The output also shows the scales for nondimensionalization and the PETSc options selected by PyLith.
This simulation did not use a fault, so PyLith used the LU preconditioner.

At the end of the output written to the termial, we see that the solver advanced the solution one time step (static simulation).
The linear solve converged after 1 iterations and the norm of the residual met the relative convergence tolerance (`ksp_rtol`) .
The nonlinear solve converged in 1 iteration, which we expect because this is a linear problem, and the residual met the absolute convergence tolerance (`snes_atol`).

## Visualizing the results

The `output` directory contains the simulation output.
Each "observer" writes its own set of files, so the solution over the domain is in one set of files, the boundary condition information is in another set of files, and the material information is in yet another set of files.
The HDF5 (`.h5`) files contain the mesh geometry and topology information along with the solution fields.
The Xdmf (`.xmf`) files contain metadata that allow visualization tools like ParaView to know where to find the information in the HDF5 files.
To visualize the data using ParaView or Visit, load the Xdmf files.

In {numref}`fig:example:reverse:2d:step01:solution` we use ParaView to visualize the displacement field using the `viz/plot_dispwarp.py` Python script.
First, we start ParaView from the `examples/reverse-2d` directory.
Next, we use the Python Shell to change the default exaggeration of the deformation to 5 to account for the large deformation.

```{code-block} python
---
caption: Change the exaggeration (warp scaling) to 5.
---
>>> WARP_SCALE = 5
```

Finally, we run the `viz/plot_dispwarp.py` Python script as described in {ref}`sec-paraview-python-scripts`.
We apply the gravitational body forces to an undeformed, stress-free domain.
As a result, the vertical deformation is about 2 kilometers.

:::{figure-md} fig:example:reverse:2d:step01:solution
<img src="figs/step01-solution.*" alt="Solution for Step 1. The colors indicate the magnitude of the displacement, and the deformation is exaggerated by a factor of 5." width="100%"/>

Solution for Step 1.
The colors of the shaded surface indicate the magnitude of the displacement, and the deformation is exaggerated by a factor of 5.
The undeformed configuration is show by the gray wireframe.
:::
