# Step 8: Use of Gravitational Body Forces

This example demonstrates the use of gravitational body forces as well as the use of initial stresses to balance the body forces. This involves enabling gravity within our domain with Dirichlet roller boundary conditions on the lateral and bottom boundaries; we do not include faults in this example.  We also demonstrate what happens when the initial stresses are not in balance with the gravitational stresses, and show how viscoelastic problems with gravitational stresses will in general not reach a steady-state solution. This example consists of three sub-problems:

* **Step 8a**: Gravitational body forces with 3-D density variations in elastic materials and initial stresses for a uniform density.
* **Step 8b**: Gravitational body forces with 3-D density variations in incompressible elastic materials.
* **Step 8c**: Gravitational body forces with 3-D density variations in elastic and viscoelastic materials and initial stresses from Step 8a plus finite strain formulation (does not reach a steady-state solution). This example is not yet implemented for PyLith v3+.

:::{warning}
Step 8c is still being updated for use with PyLith v3.
:::

:::{danger}
Step 8b takes an extremely long time to run and does not presently run in parallel.
We are working to select better preconditioners to improve the solvers and reduce the runtime.
:::

All these simulations use the default `ZeroDB` displacement boundary conditions on the lateral and bottom boundaries, and they all use `spatialdata.spatialdb.GravityField` to apply gravitational body forces. Step 8a uses reference stresses (provided in the associated spatial databases for each material) to balance the deformation induced by turning on body forces. The reference stresses assume a constant density, while there is actually a density contrast between materials, meaning that the stresses don't balance and there is some resulting deformation. Step 8b uses incompressible elasticity (using `pylith.problems.SolnDispPres`) to prevent the vertical deformation that would otherwise occur. In this case it is necessary to provide an additional Dirichlet zero pressure boundary condition on the ground surface.

{numref}`fig:example:subduction:3d:step08:diagram` shows the boundary conditions on the domain.

:::{figure-md} fig:example:subduction:3d:step08:diagram
<img src="figs/step08-diagram.*" alt="" width="100%">

Boundary conditions for gravitational body forces examples.
Body forces are applied along with roller boundary conditions on the lateral sides and bottom of the domain. For Step 8b, we also apply zero pressure BC along the upper surface of the mesh.
:::

% Features extracted from simulation parameter files.
```{include} step08a_gravity_refstate-synopsis.md
```

## Simulation parameters

The parameters specific to this example are in `step08a_gravity_refstate.cfg` and `step08b_gravity_incompressible.cfg` and include:

* `pylithapp.metadata` Metadata for this simulation. Even when the author and version are the same for all simulations in a directory, we prefer to keep that metadata in each simulation file as a reminder to keep it up-to-date for each simulation.
* `pylithapp` Parameters defining where to write the output.
* `pylithapp.problem` Parameters to set quadrature and basis order, as well as provide a gravity_field. For Step 8b, we also set the solution to `pylith.problems.SolnDispPres`.
* `pylithapp.problem.bc` For Step 8b, we add a zero pressure BC on the ground surface.
* `pylithapp.problem.materials` For Step 8a, we set `use_reference_state` and provide spatial databases that include reference stresses.

For all of the problems involving gravitational body forces we use `spatialdata.spatialdb.GravityField` to define the `gravity_field`. For Step 8a, we use higher order `quadrature_order` and `basis_order` to accurately capture the displacement field. This is not needed for Step 8b, which has much smaller displacements. Note that Step 8b also requires `mat_elastic_incompressible.cfg`, which provides all of the necessary properties using a `UniformDB`.


```{code-block} console
---
caption: Run Step 8a and 8b simulations
---
$ pylith step08a_gravity_refstate.cfg
$ pylith step08b_gravity_incompressible.cfg mat_elastic_incompressible.cfg

# The output should look something like the following. This is for Step 8b.
>> /work/charlesw/virtualenv/python3.12/lib/python3.12/site-packages/pylith/apps/PyLithApp.py:77:main
 -- pylithapp(info)
 -- Running on 1 process(es).
 >> /work/charlesw/virtualenv/python3.12/lib/python3.12/site-packages/pylith/meshio/MeshIOObj.py:38:read
 -- meshiocubit(info)
 -- Reading finite-element mesh
 >> /home/charlesw/cig/pylith3/source/pylith-fork/libsrc/pylith/meshio/MeshIOCubit.cc:148:void pylith::meshio::MeshIOCubit::_readVertices(pylith::meshio::ExodusII&, pylith::scalar_array*, int*, int*) const
 -- meshiocubit(info)
 -- Component 'reader': Reading 24824 vertices.
 >> /home/charlesw/cig/pylith3/source/pylith-fork/libsrc/pylith/meshio/MeshIOCubit.cc:208:void pylith::meshio::MeshIOCubit::_readCells(pylith::meshio::ExodusII&, pylith::int_array*, pylith::int_array*, int*, int*) const
 -- meshiocubit(info)
 -- Component 'reader': Reading 134381 cells in 4 blocks.

# -- many lines omitted --

 >> /home/charlesw/cig/pylith3/source/pylith-fork/libsrc/pylith/utils/PetscOptions.cc:239:static void pylith::utils::_PetscOptions::write(pythia::journal::info_t&, const char*, const pylith::utils::PetscOptions&)
 -- petscoptions(info)
 -- Setting PETSc options:
fieldsplit_displacement_pc_type = lu
fieldsplit_pressure_pc_type = lu
ksp_atol = 1.0e-12
ksp_converged_reason = true
ksp_error_if_not_converged = true
ksp_guess_pod_size = 8
ksp_guess_type = pod
ksp_rtol = 1.0e-12
pc_fieldsplit_schur_factorization_type = full
pc_fieldsplit_schur_precondition = full
pc_fieldsplit_type = schur
pc_type = fieldsplit
snes_atol = 1.0e-9
snes_converged_reason = true
snes_error_if_not_converged = true
snes_monitor = true
snes_rtol = 1.0e-12
ts_error_if_step_fails = true
ts_monitor = true
ts_type = beuler

 >> /work/charlesw/virtualenv/python3.12/lib/python3.12/site-packages/pylith/problems/TimeDependent.py:132:run
 -- timedependent(info)
 -- Solving problem.
0 TS dt 3.15576e+07 time 0.
    0 SNES Function norm 1.847563517818e+03
      Linear solve converged due to CONVERGED_ATOL iterations 2
    1 SNES Function norm 1.534546536971e-09
    Nonlinear solve converged due to CONVERGED_FNORM_RELATIVE iterations 1
1 TS dt 3.15576e+07 time 3.15576e+07
 >> /work/charlesw/virtualenv/python3.12/lib/python3.12/site-packages/pylith/problems/Problem.py:199:finalize
 -- timedependent(info)
 -- Finalizing problem.

```

The beginning of the output is nearly the same as in several previous examples. Both of these problems are static, however, so there is only a single time step.

## Visualizing the results

The `output` directory contains the simulation output.
Each "observer" writes its own set of files, so the solution over the domain is in one set of files, the boundary condition information is in another set of files, and the material information is in yet another set of files.
The HDF5 (`.h5`) files contain the mesh geometry and topology information along with the solution fields.
The Xdmf (`.xmf`) files contain metadata that allow visualization tools like ParaView to know where to find the information in the HDF5 files.
To visualize the data using ParaView or Visit, load the Xdmf files.

In {numref}`fig:example:subduction:3d:step08a:solution` we use ParaView to visualize the x displacement field using the `viz/plot_dispwarp.py` Python script.
We start ParaView from the `examples/subduction-3d` directory and then run the `viz/plot_dispwarp.py` Python script as described in {ref}`sec-paraview-python-scripts`.

:::{figure-md} fig:example:subduction:3d:step08a:solution
<img src="figs/step08a-solution.*" alt="Solution for Step 8a. The colors indicate the z displacement, and the deformation is exaggerated by a factor of 50." width="100%"/>

Solution for Step 8a.
The colors of the shaded surface indicate the x displacement, and the deformation is exaggerated by a factor of 50.
:::

Note that the vertical displacements for Step 8b are very close to zero, because the material is nearly incompressible.
