# Step 1: Static Coseismic Slip

% Features extracted from simulation parameter files.
```{include} step01_coseismic-synopsis.md
```

## Simulation parameters

This example involves a static simulation that solves for the deformation from prescribed coseismic slip on the subduction interface.
The depth variation in the prescribed slip is based on the 2011 Tohoku-oki earthquake.
{numref}`fig:example:subduction:2d:step01:diagram` shows the boundary conditions on the domain.
The parameters specific to this example are in `step01_coseismic.cfg`.

:::{figure-md} fig:example:subduction:2d:step01:diagram
<img src="figs/step01-diagram.*" alt="" width="100%">

Boundary conditions for static coseismic slip on the subduction interface.
We prescribe reverse slip that varies with depth and roller boundary conditions on the lateral sides and bottom of the domain.
:::

We only presceibe slip on the subduction interface, so we create an array with 1 fault.
We specify slip as a function of depth, so we use a `SimpleDB` with linear interpolation.

```{code-block} cfg
---
caption: Fault parameters for Step 1.
---
[pylithapp.problem]
interfaces = [fault]

[pylithapp.problem.interfaces.fault]
label = fault_slabtop
label_value = 21
edge = fault_slabtop_edge
edge_value = 31

observers.observer.data_fields = [slip]

[pylithapp.problem.interfaces.fault.eq_ruptures.rupture]
db_auxiliary_field = spatialdata.spatialdb.SimpleDB
db_auxiliary_field.description = Fault rupture auxiliary field spatial database
db_auxiliary_field.iohandler.filename = fault_coseismic.spatialdb
db_auxiliary_field.query_type = linear
```

```{code-block} cfg
---
caption: Dirichlet boundary condition parameters for Step 1. We only show the details for the east east boundary of the crust.
---
[pylithapp.problem]
bc = [bc_east_crust, bc_east_mantle, bc_west, bc_bottom]

[pylithapp.problem.bc.bc_east_crust]
label = bndry_east_crust
label_value = 12
constrained_dof = [0]
db_auxiliary_field = pylith.bc.ZeroDB
db_auxiliary_field.description = Dirichlet BC on east boundary (crust)
```

## Running the simulation

```{code-block} console
---
caption: Run Step 1 simulation
---
$ pylith step01_coseismic.cfg

# The output should look something like the following.
 >> /software/unix/py39-venv/pylith-debug/lib/python3.9/site-packages/pylith/meshio/MeshIOObj.py:44:read
 -- meshiopetsc(info)
 -- Reading finite-element mesh
 >> /src/cig/pylith/libsrc/pylith/meshio/MeshIO.cc:94:void pylith::meshio::MeshIO::read(topology::Mesh *)
 -- meshiopetsc(info)
 -- Component 'reader': Domain bounding box:
    (-600000, 600000)
    (-600000, 399.651)
 >> /software/unix/py39-venv/pylith-debug/lib/python3.9/site-packages/pylith/problems/Problem.py:116:preinitialize
 -- timedependent(info)
 -- Performing minimal initialization before verifying configuration.
 >> /software/unix/py39-venv/pylith-debug/lib/python3.9/site-packages/pylith/problems/Solution.py:44:preinitialize
 -- solution(info)
 -- Performing minimal initialization of solution.
 >> /software/unix/py39-venv/pylith-debug/lib/python3.9/site-packages/pylith/problems/Problem.py:175:verifyConfiguration
 -- timedependent(info)
 -- Verifying compatibility of problem configuration.
 >> /software/unix/py39-venv/pylith-debug/lib/python3.9/site-packages/pylith/problems/Problem.py:221:_printInfo
 -- timedependent(info)
 -- Scales for nondimensionalization:
    Length scale: 1000*m
    Time scale: 3.15576e+09*s
    Pressure scale: 3e+10*m**-1*kg*s**-2
    Density scale: 2.98765e+23*m**-3*kg
    Temperature scale: 1*K
 >> /software/unix/py39-venv/pylith-debug/lib/python3.9/site-packages/pylith/problems/Problem.py:186:initialize
 -- timedependent(info)
 -- Initializing timedependent problem with quasistatic formulation.
 >> /src/cig/pylith/libsrc/pylith/utils/PetscOptions.cc:235:static void pylith::utils::_PetscOptions::write(pythia::journal::info_t &, const char *, const pylith::utils::PetscOptions &)
 -- petscoptions(info)
 -- Setting PETSc options:
fieldsplit_displacement_ksp_type = preonly
fieldsplit_displacement_pc_type = lu
fieldsplit_lagrange_multiplier_fault_ksp_type = preonly
fieldsplit_lagrange_multiplier_fault_pc_type = lu
ksp_atol = 1.0e-12
ksp_converged_reason = true
ksp_error_if_not_converged = true
ksp_rtol = 1.0e-12
pc_fieldsplit_schur_factorization_type = lower
pc_fieldsplit_schur_precondition = selfp
pc_fieldsplit_schur_scale = 1.0
pc_fieldsplit_type = schur
pc_type = fieldsplit
pc_use_amat = true
snes_atol = 1.0e-9
snes_converged_reason = true
snes_error_if_not_converged = true
snes_monitor = true
snes_rtol = 1.0e-12
ts_error_if_step_fails = true
ts_monitor = true
ts_type = beuler

 >> /software/unix/py39-venv/pylith-debug/lib/python3.9/site-packages/pylith/problems/TimeDependent.py:139:run
 -- timedependent(info)
 -- Solving problem.
0 TS dt 0.05 time -0.05
    0 SNES Function norm 5.454651006059e-01 
    Linear solve converged due to CONVERGED_ATOL iterations 85
    1 SNES Function norm 3.437540896787e-12 
  Nonlinear solve converged due to CONVERGED_FNORM_ABS iterations 1
1 TS dt 0.05 time 0.
 >> /software/unix/py39-venv/pylith-debug/lib/python3.9/site-packages/pylith/problems/Problem.py:201:finalize
 -- timedependent(info)
 -- Finalizing problem.
```

At the beginning of the output written to the terminal, we see that PyLith is reading the mesh using the `MeshIOPetsc` reader and that it found the domain to extend from -600 km to +600 km in the x direction and from -600 km to 0 in the y direction.
The output also includes the scales used for nondimensionalization and the default PETSc options.

At the end of the output written to the terminal, we see that the solver advanced the solution one time step (static simulation).
The linear solve converged after 85 iterations and the norm of the residual met the absolute convergence tolerance (`ksp_atol`) .
The nonlinear solve converged in 1 iteration, which we expect because this is a linear problem, and the residual met the absolute convergence tolerance (`snes_atol`).

## Visualizing the results

The `output` directory contains the simulation output.
Each "observer" writes its own set of files, so the solution over the domain is in one set of files, the boundary condition information is in another set of files, and the material information is in yet another set of files.
The HDF5 (`.h5`) files contain the mesh geometry and topology information along with the solution fields.
The Xdmf (`.xmf`) files contain metadata that allow visualization tools like ParaView to know where to find the information in the HDF5 files.
To visualize the data using ParaView or Visit, load the Xdmf files.

In {numref}`fig:example:subduction:2d:step01:solution` we use ParaView to visualize the x displacement field using the `viz/plot_dispwarp.py` Python script.
We start ParaView from the `examples/subduction-2d` directory and then run the `viz/plot_dispwarp.py` Python script as described in {ref}`sec-paraview-python-scripts`.

:::{figure-md} fig:example:subduction:2d:step01:solution
<img src="figs/step01-solution.*" alt="Solution for Step 1. The colors indicate the magnitude of the x displacement, and the deformation is exaggerated by a factor of 1000." width="100%"/>

Solution for Step 1.
The colors of the shaded surface indicate the magnitude of the x displacement, and the deformation is exaggerated by a factor of 1000.
:::
