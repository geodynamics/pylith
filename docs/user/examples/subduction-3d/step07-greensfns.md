# Step 7: Generation of Green's Functions and Slow Slip Inversion

In this example we compute the static Green's function for prescribed slip impulses on the central portion of the subduction interface.
We then use the responses to perform a simple geodetic inversion of the slip from Step 6.
We divide generating Green's functions into two sub-problems:

* **Step 7a**: Left-lateral slip component
* **Step 7b**: Reverse slip component

Although PyLith can generate the two components in one simulation, we often prefer to speed up the process by running simulations for each of the components simultaneously using multiple processes on a cluster.

{numref}`fig:example:subduction:3d:step07:diagram` shows the boundary conditions on the domain.

:::{figure-md} fig:example:subduction:3d:step07:diagram
<img src="figs/step07-diagram.*" alt="" width="100%">

Boundary conditions for Green's function generation.
We apply left-lateral (Step 7a) and reverse (Step 7b) slip impulses on a patch within the subduction interface with roller boundary conditions on the lateral sides and bottom of the domain.
:::

% Features extracted from simulation parameter files.
```{include} step07a_leftlateral-synopsis.md
```

## Simulation parameters

The parameters specific to this example are in `step07a_leftlateral.cfg` and `step07b_reverse.cfg` and include:

* `pylithapp.metadata` Metadata for this simulation. Even when the author and version are the same for all simulations in a directory, we prefer to keep that metadata in each simulation file as a reminder to keep it up-to-date for each simulation.
* `pylithapp` Parameters defining where to write the output. We also change the problem type to `pylith.problems.GreensFns`.
* `pylithapp.problem` Parameters for the time step information as well as solution field with displacement and Lagrange multiplier subfields.
* `pylithapp.journal.info` We turn on `greensfns` info for this problem.
* `pylithapp.interfaces` Parameters for the slip impulses applied on the top of the slab.

For Green's functions we use the `FaultCohesiveImpulses` kinematic source to prescribe slip impulses for left-lateral (Step 7a -- `impulse_dof = [1]`) and reverse (Step 7b -- `impulse_dof = [2]`) slip directions. As for Step 6, we use `OutputSolnPoints` to output the solution at fake CGNSS stations.

```{code-block} console
---
caption: Run Step 7a and 7b simulations
---
$ pylith step07a_leftlateral.cfg mat_elastic.cfg
$ pylith step07b_reverse.cfg mat_elastic.cfg

# The output should look something like the following.
>> software/virtualenv/python312/lib/python3.12/site-packages/pylith/apps/PyLithApp.py:77:main
 -- pylithapp(info)
 -- Running on 1 process(es).
 >> software/virtualenv/python312/lib/python3.12/site-packages/pylith/meshio/MeshIOObj.py:38:read
 -- meshiocubit(info)
 -- Reading finite-element mesh
 >> software/cig/pylith3/source/pylith-fork/libsrc/pylith/meshio/MeshIOCubit.cc:148:void pylith::meshio::MeshIOCubit::_readVertices(ExodusII &, scalar_array *, int *, int *) const
 -- meshiocubit(info)
 -- Component 'reader': Reading 24738 vertices.
 >> software/cig/pylith3/source/pylith-fork/libsrc/pylith/meshio/MeshIOCubit.cc:208:void pylith::meshio::MeshIOCubit::_readCells(ExodusII &, int_array *, int_array *, int *, int *) const
 -- meshiocubit(info)
 -- Component 'reader': Reading 133827 cells in 4 blocks.

# -- many lines omitted --

>> software/cig/pylith3/source/pylith-fork/libsrc/pylith/utils/PetscOptions.cc:239:static void pylith::utils::_PetscOptions::write(pythia::journal::info_t &, const char *, const PetscOptions &)
 -- petscoptions(info)
 -- Setting PETSc options:
dm_reorder_section = true
dm_reorder_section_type = cohesive
ksp_atol = 1.0e-12
ksp_converged_reason = true
ksp_error_if_not_converged = true
ksp_guess_pod_size = 8
ksp_guess_type = pod
ksp_rtol = 1.0e-12
mg_fine_pc_type = vpbjacobi
pc_type = gamg
snes_atol = 1.0e-9
snes_converged_reason = true
snes_error_if_not_converged = true
snes_monitor = true
snes_rtol = 1.0e-12
ts_error_if_step_fails = true
ts_monitor = true
ts_type = beuler

# -- many lines omitted --

>> software/cig/pylith3/source/pylith-fork/libsrc/pylith/problems/GreensFns.cc:322:void pylith::problems::GreensFns::solve()
 -- greensfns(info)
 -- Component 'problem': Computing Green's function 92 of 93.
  0 SNES Function norm 7.551182042932e-02
    Linear solve converged due to CONVERGED_ATOL iterations 20
  1 SNES Function norm 1.971805442527e-10
  Nonlinear solve converged due to CONVERGED_ITS iterations 1
 >> software/cig/pylith3/source/pylith-fork/libsrc/pylith/problems/GreensFns.cc:322:void pylith::problems::GreensFns::solve()
 -- greensfns(info)
 -- Component 'problem': Computing Green's function 93 of 93.
  0 SNES Function norm 6.572171502838e-02
    Linear solve converged due to CONVERGED_ATOL iterations 21
  1 SNES Function norm 1.378742924611e-10
  Nonlinear solve converged due to CONVERGED_ITS iterations 1
 >> software/virtualenv/python312/lib/python3.12/site-packages/pylith/problems/Problem.py:199:finalize
 -- greensfns(info)
 -- Finalizing problem.
```

The beginning of the output is nearly the same as in several previous examples.
In this case, however, rather than stepping through time, we are stepping through impulse numbers. In each simulation, there are 93 impulses applied on the fault.

## Visualizing the results

The `output` directory contains the simulation output.
Each "observer" writes its own set of files, so the solution over the domain is in one set of files, the boundary condition information is in another set of files, and the material information is in yet another set of files.
The HDF5 (`.h5`) files contain the mesh geometry and topology information along with the solution fields.
The Xdmf (`.xmf`) files contain metadata that allow visualization tools like ParaView to know where to find the information in the HDF5 files.
To visualize the data using ParaView or Visit, load the Xdmf files.

In {numref}`fig:example:subduction:3d:step07a:slip` we use ParaView to visualize the x displacement field using the `viz/plot_dispwarp.py` Python script.
We start ParaView from the `examples/subduction-3d` directory and then run the `viz/plot_dispwarp.py` Python script as described in {ref}`sec-paraview-python-scripts`.

:::{figure-md} fig:example:subduction:3d:step07a:slip
<img src="figs/step07a-slipimpulse.*" alt="Slip impulse 41 for Step 7a. The colors indicate the magnitude of the slip." width="100%"/>

Distribution of left-lateral fault slip for slip impulse 41 (time corresponds to the zero-based index) for Step 7a.
The colors of the shaded surface indicate the magnitude of the slip.
:::

## Performing a simple inversion

To perform our simulated inversion, we first postprocess the displacement results from Step 6.
We do this using the `make_synthetic_gnssdisp.py` Python script, which in turn reads the parameters from the `make_synthetic_gnssdisp.cfg` file:

```{code-block} console
---
caption: Generate synthetic data for inversion
---
$ ./make_synthetic_gnssdisp.py
```
This script will produce the files `cgnss_synthetic_displacement.txt` and `cgnss_synthetic_displacement.vtk`, which are used in the inversion.
Once these have been created (and both sets of Green's functions have been generated), we can perform a simple inversion using the `slip_invert.py` script, which in turn reads the inversion parameters from `slip_invert.cfg`.

```{code-block} console
---
caption: Perform inversion of synthetic data
---
$ ./slip_invert.py
```

The Python script performs a simple weighted generalized least-squares inversion, with penalty values applied to the deviation of the solution from the a priori value (zero, in this case).
The slip inversion will produce several files in the `output` directory:

* `step07-inversion-summary.txt` A summary of the model residuals for each penalty parameter.
* `step07-inversion-displacement.h5` The predicted site displacements for different penalty parameters.
* `step07-inversion-slip.h5` The predicted fault slip for different penalty parameters.

Note that the inversion produces a `.xmf` file for each of the HDF5 files.
We provide a Matplotlib Python script to visualize the data misfit for various penalty parameters:

```{code-block} console
---
caption: Visualize data misfit
---
$ viz/plot_inversion_misfit.py --summary=output/step07-inversion-summary.txt
```

{numref}`fig:example:subduction:3d:step07:inversion_summary` shows the data residual versus the penalty residual, showing very clearly the 'corner' of the L-curve.

:::{figure-md} fig:example:subduction:3d:step07:inversion_summary
<img src="figs/step07-inversion-summary.*" alt="" width="75%">

Data residual versus penalty residual for inversion from Step 7.
:::