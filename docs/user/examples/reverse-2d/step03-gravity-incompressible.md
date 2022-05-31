# Step 3: Gravitational Body Forces with Incompressible Elasticity

In this example we use incompressible elasticity (see {ref}`sec-user-governing-eqns-incompressible-elasticity` for the governing equations and finite-element formulation) to obtain the stress field associated with gravitational body forces.
Because the material is incompressible and the material is confined on the lateral boundaries and bottom, we do not expect any deformation.
In general, this is a more robust way to determine an initial stress state for gravitational body forces compared to using a reference stress state, especially when the material properties are not uniform.
We use the same roller boundary conditions that we used in Steps 1 and 2.

% Metadata extracted from parameter files.
```{include} step03_gravity_incompressible-synopsis.md
```

## Simulation parameters

The parameters specific to this example are in `step03_gravity_incompressible.cfg`.
These include:

* `pylithapp.metadata` Metadata for this simulation. Even when the author and version are the same for all simulations in a directory, we prefer to keep that metadata in each simulation file as a reminder to keep it up-to-date for each simulation.
* `pylithapp` Parameters defining where to write the output.
* `pylithapp.problem` Parameters for specifying the gravitational body forces and using displacement and pressure in the solution field for incompressible elasticity.
* `pylithapp.problem.materials` Use linear, isotropic incompressible elasticity.
* `pylithapp.problem.bc` Add a Dirichlet boundary condition to specify zero pressure on the top surface (+y boundary).

```{code-block} console
---
caption: Run Step 3 simulation
---
$ pylith step03_gravity_incompressible.cfg

# The output should look something like the following.
 >> /software/baagaard/py38-venv/pylith-debug/lib/python3.8/site-packages/pylith/meshio/MeshIOObj.py:44:read
 -- meshiopetsc(info)
 -- Reading finite-element mesh
 >> /home/baagaard/src/cig/pylith/libsrc/pylith/meshio/MeshIO.cc:94:void pylith::meshio::MeshIO::read(pylith::topology::Mesh*)
 -- meshiopetsc(info)
 -- Component 'reader': Domain bounding box:
    (-100000, 100000)
    (-100000, 0)

# -- many lines omitted --

 >> /Users/baagaard/src/cig/pylith/libsrc/pylith/utils/PetscOptions.cc:235:static void pylith::utils::_PetscOptions::write(pythia::journal::info_t &, const char *, const pylith::utils::PetscOptions &)
 -- petscoptions(info)
 -- Setting PETSc options:
ksp_atol = 1.0e-12
ksp_converged_reason = true
ksp_error_if_not_converged = true
ksp_rtol = 1.0e-12
pc_fieldsplit_schur_factorization_type = full
snes_atol = 1.0e-9
snes_converged_reason = true
snes_error_if_not_converged = true
snes_monitor = true
snes_rtol = 1.0e-12
ts_error_if_step_fails = true
ts_monitor = true
ts_type = beuler

 >> /Users/baagaard/src/cig/pylith/libsrc/pylith/utils/PetscOptions.cc:235:static void pylith::utils::_PetscOptions::write(pythia::journal::info_t &, const char *, const pylith::utils::PetscOptions &)
 -- petscoptions(info)
 -- Ignoring PETSc options (already set):
fieldsplit_displacement_pc_type = lu
fieldsplit_pressure_pc_type = lu
pc_fieldsplit_schur_precondition = full
pc_fieldsplit_type = schur
pc_type = fieldsplit

 >> /Users/baagaard/software/unix/py39-venv/pylith-debug/lib/python3.9/site-packages/pylith/problems/TimeDependent.py:139:run
 -- timedependent(info)
 -- Solving problem.
0 TS dt 0.01 time 0.
    0 SNES Function norm 4.866941773461e-01 
    Linear solve converged due to CONVERGED_ATOL iterations 1
    1 SNES Function norm 3.099989574301e-13 
  Nonlinear solve converged due to CONVERGED_FNORM_ABS iterations 1
1 TS dt 0.01 time 0.01
 >> /Users/baagaard/software/unix/py39-venv/pylith-debug/lib/python3.9/site-packages/pylith/problems/Problem.py:201:finalize
 -- timedependent(info)
 -- Finalizing problem.
```

PyLith detected use of incompressible elasticity, so it selected a field split preconditioner with an LU preconditioner for each of the solution subfields as described in {ref}`sec-user-run-pylith-petsc-options`.
As a result, the linear solve converged in 1 iterations.

## Visualizing the results

The `output` directory contains the simulation output.
Each "observer" writes its own set of files, so the solution over the domain is in one set of files, the boundary condition information is in another set of files, and the material information is in yet another set of files.
The HDF5 (`.h5`) files contain the mesh geometry and topology information along with the solution fields.
The Xdmf (`.xmf`) files contain metadata that allow visualization tools like ParaView to know where to find the information in the HDF5 files.
To visualize the data using ParaView or Visit, load the Xdmf files.

In {numref}`fig:example:reverse:2d:step03:solution` we use ParaView to visualize the displacement field using the `viz/plot_dispwarp.py` Python script.
First, we start ParaView from the `examples/reverse-2d` directory.
Before running the `viz/plot_dispwarp.py` Python script as described in {ref}`sec-paraview-python-scripts`, we set the simulation name in the ParaView Python Shell.

```{code-block} python
---
caption: Set the simulation in the ParaView Python Shell.
---
>>> SIM = "step03_gravity_incompressible"
```

:::{figure-md} fig:example:reverse:2d:step03:solution
<img src="figs/step03-solution.*" alt="Solution for Step 3. The colors indicate the magnitude of the displacement." width="100%"/>

Solution for Step 3.
The colors of the shaded surface indicate the magnitude of the displacement.
The undeformed configuration is show by the gray wireframe.
There is negligible deformation and the stress state (not shown) matches the one in Step 2.
:::
