# Step 5: Static Coseismic Slip

This example involves a static simulation that solves for the deformation from prescribed coseismic slip on the main fault.
We specify 2 meters of reverse slip.
{numref}`fig:example:reverse:2d:step05:diagram` shows the boundary conditions on the domain.

:::{figure-md} fig:example:reverse:2d:step05:diagram
<img src="figs/step05-diagram.*" alt="" scale="75%">

Boundary conditions for static coseismic slip on the main fault.
We prescribe 2 meters of reverse slip with roller boundary conditions on the lateral sides and bottom of the domain.
:::

% Metadata extracted from parameter files.
```{include} step05_onefault-synopsis.md
```

## Simulation parameters

The parameters specific to this example are in `step05_onefault.cfg` and include:

* `pylithapp.metadata` Metadata for this simulation. Even when the author and version are the same for all simulations in a directory, we prefer to keep that metadata in each simulation file as a reminder to keep it up-to-date for each simulation.
* `pylithapp` Parameters defining where to write the output.
* `pylithapp.problem` Parameters for the solution field with displacement and Lagrange multiplier subfields.
* `pylithapp.problem.fault` Parameters for prescribed slip on the fault.

:::{important}
In 2D simulations slip is specified in terms of opening and left-lateral components.
This provides a consistent, unique sense of slip that is independent of the fault orientation.
For our geometry in this example, right lateral slip corresponds to reverse slip on the dipping fault.
:::

:::{important}
The main fault contains one end that is buried within the domain.
When PyLith inserts cohesive cells into a mesh with buried edges (in this case a point), we must identify these buried edges so that PyLith properly adjusts the topology along these edges.
:::

```{code-block} console
---
caption: Run Step 5 simulation
---
$ pylith step05_onefault.cfg

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

 >> /Users/baagaard/software/unix/py39-venv/pylith-debug/lib/python3.9/site-packages/pylith/problems/TimeDependent.py:139:run
 -- timedependent(info)
 -- Solving problem.
0 TS dt 0.01 time 0.
    0 SNES Function norm 2.129295960330e-02 
    Linear solve converged due to CONVERGED_ATOL iterations 330
    1 SNES Function norm 6.693855246577e-13 
  Nonlinear solve converged due to CONVERGED_FNORM_ABS iterations 1
1 TS dt 0.01 time 0.01
 >> /Users/baagaard/software/unix/py39-venv/pylith-debug/lib/python3.9/site-packages/pylith/problems/Problem.py:201:finalize
 -- timedependent(info)
 -- Finalizing problem.
```

From the end of the output written to the terminal window, we see that the linear solver converged in 330 iterations and met the absolute convergence tolerance (`ksp_atol`).
As we expect for this linear problem, the nonlinear solver converged in 1 iteration.
You can run this problem using the default PETSc options for a parallel simulation and see how the use of the algebraic multigrid preconditioner reduces the number of iterations.

```{code-block} console
---
caption: Run Step 5 simulation
---
$ pylith step05_onefault.cfg --problem.petsc_defaults.parallel

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
fieldsplit_displacement_ksp_type = preonly
fieldsplit_displacement_mg_levels_ksp_type = richardson
fieldsplit_displacement_mg_levels_pc_type = sor
fieldsplit_displacement_pc_type = gamg
fieldsplit_lagrange_multiplier_fault_ksp_type = preonly
fieldsplit_lagrange_multiplier_fault_mg_levels_ksp_type = richardson
fieldsplit_lagrange_multiplier_fault_mg_levels_pc_type = sor
fieldsplit_lagrange_multiplier_fault_pc_type = gamg
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

 >> /Users/baagaard/software/unix/py39-venv/pylith-debug/lib/python3.9/site-packages/pylith/problems/TimeDependent.py:139:run
 -- timedependent(info)
 -- Solving problem.
0 TS dt 0.01 time 0.
    0 SNES Function norm 2.129295960330e-02 
    Linear solve converged due to CONVERGED_ATOL iterations 47
    1 SNES Function norm 1.185726544112e-12 
  Nonlinear solve converged due to CONVERGED_FNORM_ABS iterations 1
1 TS dt 0.01 time 0.01
 >> /Users/baagaard/software/unix/py39-venv/pylith-debug/lib/python3.9/site-packages/pylith/problems/Problem.py:201:finalize
 -- timedependent(info)
 -- Finalizing problem.
```

Notice the use of the `GAMG` preconditioner for the displacement and Lagrange multiplier solution subfields.

## Visualizing the results

The `output` directory contains the simulation output.
Each "observer" writes its own set of files, so the solution over the domain is in one set of files, the boundary condition information is in another set of files, and the material information is in yet another set of files.
The HDF5 (`.h5`) files contain the mesh geometry and topology information along with the solution fields.
The Xdmf (`.xmf`) files contain metadata that allow visualization tools like ParaView to know where to find the information in the HDF5 files.
To visualize the data using ParaView or Visit, load the Xdmf files.

In {numref}`fig:example:reverse:2d:step05:solution` we use ParaView to visualize the x displacement field using the `viz/plot_dispwarp.py` Python script.
First, we start ParaView from the `examples/reverse-2d` directory.
Before running the `viz/plot_dispwarp.py` Python script as described in {ref}`sec-paraview-python-scripts`, we set the simulation name in the ParaView Python Shell.

```{code-block} python
---
caption: Set the simulation and exaggeration in the ParaView Python Shell.
---
>>> SIM = "step05_onefault"
>>> FIELD_COMPONENT = "X"
```

:::{figure-md} fig:example:reverse:2d:step05:solution
<img src="figs/step05-solution.*" alt="Solution for Step 5. The colors indicate the magnitude of the x displacement, and the deformation is exaggerated by a factor of 1000." width="100%"/>

Solution for Step 5.
The colors of the shaded surface indicate the magnitude of the x displacement, and the deformation is exaggerated by a factor of 1000.
The undeformed configuration is show by the gray wireframe.
:::
