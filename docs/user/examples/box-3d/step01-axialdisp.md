# Step 1: Axial Extension

% Meatadata extracted from parameter files
```{include} step01_axialdisp-synopsis.md
```

## Simulation parameters

This example corresponds to axial extension in the x direction and axial compression in the y direction.
We apply Dirichlet (displacement) boundary conditions for the x displacement on the +x (`boundary_xpos`) and -x (`boundary_xneg`) boundaries and for the y displacement on the +y (`boundary_ypos`) and -y (`boundary_yneg`) boundaries.
We apply roller Dirichlet boundary conditions on the -z (`boundary_zneg`) boundary.
{numref}`fig:example:box:3d:step01:diagram` shows the boundary conditions on the domain.
The parameters specific to this example are in `step01_axialdisp.cfg`.

:::{figure-md} fig:example:box:3d:step01:diagram
<img src="figs/step01-diagram.*" alt="" scale="75%">

Boundary conditions for axial extension in the x direction and axial compression in the y direction.
We constrain the x displacement on the +x and -x boundaries, the y displacement on the +y and -y boundaries, and set the z displacement to zero on the -z boundary.
:::

We create an array of 5 `DirichletTimeDependent` boundary conditions.
For each of these boundary conditions we must specify which degrees of freedom are constrained, the name of the label marking the boundary (name of the group of vertices in the finite-element mesh file), and the values for the Dirichlet boundary condition.

```{code-block} cfg
---
caption: Specifying the boundary conditions for Step 1. We only show the detailed settings for the +x boundary.
---
[pylithapp.problem]
bc = [bc_xneg, bc_xpos, bc_yneg, bc_ypos, bc_zneg]
bc.bc_xneg = pylith.bc.DirichletTimeDependent
bc.bc_xpos = pylith.bc.DirichletTimeDependent
bc.bc_yneg = pylith.bc.DirichletTimeDependent
bc.bc_ypos = pylith.bc.DirichletTimeDependent
bc.bc_zneg = pylith.bc.DirichletTimeDependent
        
[pylithapp.problem.bc.bc_xpos]
label = boundary_xpos
label_value = 11
constrained_dof = [0]

db_auxiliary_field = spatialdata.spatialdb.UniformDB
db_auxiliary_field.description = Dirichlet BC on +x boundary
db_auxiliary_field.values = [initial_amplitude_x, initial_amplitude_y, initial_amplitude_z]
db_auxiliary_field.data = [+2.0*m, 0*m, 0*m]
```

## Running the simulation

```{code-block} console
---
caption: Run Step 1 simulation
---
$ pylith step01_axialdisp.cfg

# The output should look something like the following.
 >> /software/unix/py39-venv/pylith-debug/lib/python3.9/site-packages/pylith/meshio/MeshIOObj.py:44:read
 -- meshiopetsc(info)
 -- Reading finite-element mesh
 >> /src/cig/pylith/libsrc/pylith/meshio/MeshIO.cc:94:void pylith::meshio::MeshIO::read(topology::Mesh *)
 -- meshiopetsc(info)
 -- Component 'reader': Domain bounding box:
    (-6000, 6000)
    (-6000, 6000)
    (-9000, 0)
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
    Length scale: 100000*m
    Displacement scale: 1*m
    Time scale: 3.15576e+09*s
    Rigidity scale: 1e+10*m**-1*kg*s**-2
    Temperature scale: 1*K
 >> /software/unix/py39-venv/pylith-debug/lib/python3.9/site-packages/pylith/problems/Problem.py:186:initialize
 -- timedependent(info)
 -- Initializing timedependent problem with quasistatic formulation.
 >> /src/cig/pylith/libsrc/pylith/utils/PetscOptions.cc:235:static void pylith::utils::_PetscOptions::write(pythia::journal::info_t &, const char *, const pylith::utils::PetscOptions &)
 -- petscoptions(info)
 -- Setting PETSc options:
ksp_atol = 1.0e-7
ksp_converged_reason = true
ksp_error_if_not_converged = true
ksp_guess_pod_size = 8
ksp_guess_type = pod
ksp_rtol = 1.0e-12
mg_fine_ksp_max_it = 5
mg_levels_pc_type = pbjacobi
pc_gamg_coarse_eq_limit = 200
pc_type = gamg
snes_atol = 4.0e-7
snes_converged_reason = true
snes_error_if_not_converged = true
snes_monitor = true
snes_rtol = 1.0e-12
ts_error_if_step_fails = true
ts_exact_final_time = matchstep
ts_monitor = true
ts_type = beuler
viewer_hdf5_collective = true

 >> /software/unix/py39-venv/pylith-debug/lib/python3.9/site-packages/pylith/problems/TimeDependent.py:139:run
 -- timedependent(info)
 -- Solving problem.
0 TS dt 0.001 time 0.
    0 SNES Function norm 2.045040813383e+00
      Linear solve converged due to CONVERGED_ATOL iterations 4
    1 SNES Function norm 3.894554358059e-09
    Nonlinear solve converged due to CONVERGED_FNORM_ABS iterations 1
1 TS dt 0.001 time 0.001
 >> /software/unix/py39-venv/pylith-debug/lib/python3.9/site-packages/pylith/problems/Problem.py:201:finalize
 -- timedependent(info)
 -- Finalizing problem.
WARNING! There are options you set that were not used!
WARNING! could be spelling mistake, etc!
There is one unused database option. It is:
Option left: name:-mg_levels_pc_type value: pbjacobi source: code
```

At the beginning of the output written to the terminal, we see that PyLith is reading the mesh using the `MeshIOAscii` reader and that it found the domain to extend from -6000 m to +6000 m in the x and y directions and from -9000 m to 0 m in the z direction.
We also see the scales used to nondimensionalize the problem.

Near the end of the output written to the terminal, we see the PETSc options PyLith selected based on the governing equations and formulation as discussed in {ref}`sec-user-run-pylith-petsc-options`.
The solver advanced the solution one time step (static simulation).
The linear solve converged in 4 iterations.
The norm of the residual met the absolute tolerance convergence criterion (`ksp_atol`).
The nonlinear solve converged in 1 iteration, which we expect because this is a linear problem, and the residual met the absolute convergence tolerance (`snes_atol`).
For this small problem, the multi-grid preconditioner has fewer levels, so we get a warning about an unused PETSc option.

## Visualizing the results

The `output` directory contains the simulation output.
Each "observer" writes its own set of files, so the solution over the domain is in one set of files, the boundary condition information is in another set of files, and the material information is in yet another set of files.
The HDF5 (`.h5`) files contain the mesh geometry and topology information along with the solution fields.
The Xdmf (`.xmf`) files contain metadata that allow visualization tools like ParaView to know where to find the information in the HDF5 files.
To visualize the data using ParaView or Visit, load the Xdmf files.

In {numref}`fig:example:box:3d:step01:solution` we use the `pylith_viz` utility to visualize the x displacement field.

```{code-block} console
---
caption: Visualize PyLith output using `pylith_viz`.
---
pylith_viz --filenames=output/step01_axialdisp-domain.h5 warp_grid --component=x
```

:::{figure-md} fig:example:box:3d:step01:solution
<img src="figs/step01-solution.*" alt="Solution for Step 1. The colors indicate the x displacement, and the deformation is exaggerated by a factor of 1000." width="600px"/>

Solution for Step 1.
The colors of the shaded surface indicate the x displacement, and the deformation is exaggerated by a factor of 1000.
The undeformed configuration is shown by the gray wireframe.
:::
