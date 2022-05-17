# Step 2: Gravitational Body Forces with Reference Stress

This example involves using a reference stress state to minimize the deformation when we apply the gravitational body forces.
The solution will be the perturbation from the reference state with zero displacements.
This is one method for obtaining an initial stress state associated with gravitational body forces.
We use the same roller boundary conditions that we used in Step 1.

% Metadata extracted from parameter files.
```{include} step02_gravity_refstate-synopsis.md
```

## Simulation parameters

The parameters specific to this example are in `step02_gravity_refstate.cfg`.
These include:

* `pylithapp.metadata` Metadata for this simulation. Even when the author and version are the same for all simulations in a directory, we prefer to keep that metadata in each simulation file as a reminder to keep it up-to-date for each simulation.
* `pylithapp` Parameters defining where to write the output.
* `pylithapp.problem` Parameters for specifying the gravitational body forces and adjusting the basis order.
* `pylithapp.problem.materials` Parameters for setting the reference stress state.

We use a reference stress state that matches the overburden (lithostatic) pressure.
We have uniform material properties, so the overbudern is
%
\begin{equation}
\sigma_{xx} = \sigma_{yy} = \sigma_{zz} = \int_0^z \rho g \, dz = \rho g z,
\end{equation}
where compressive stress is negative.

```{code-block} console
---
caption: Run Step 2 simulation
---
$ pylith step02_gravity_refstate.cfg

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
    0 SNES Function norm 4.578015693966e-15 
  Nonlinear solve converged due to CONVERGED_FNORM_ABS iterations 0
1 TS dt 0.01 time 0.01
 >> /Users/baagaard/software/unix/py39-venv/pylith-debug/lib/python3.9/site-packages/pylith/problems/Problem.py:201:finalize
 -- timedependent(info)
 -- Finalizing problem.
WARNING! There are options you set that were not used!
WARNING! could be spelling mistake, etc!
There is one unused database option. It is:
Option left: name:-ksp_converged_reason (no value)
```

By design we set the reference stress state so that it matches the loading from gravitational body forces in our domain with uniform material properties.
As a result, the first nonlinear solver residual evaluation meets the convergence criteria.
The linear solver is not used; this is why PETSc reports an unused option at the end of the simulation.

## Visualizing the results

The `output` directory contains the simulation output.
Each "observer" writes its own set of files, so the solution over the domain is in one set of files, the boundary condition information is in another set of files, and the material information is in yet another set of files.
The HDF5 (`.h5`) files contain the mesh geometry and topology information along with the solution fields.
The Xdmf (`.xmf`) files contain metadata that allow visualization tools like ParaView to know where to find the information in the HDF5 files.
To visualize the data using ParaView or Visit, load the Xdmf files.

In {numref}`fig:example:reverse:2d:step02:solution` we use ParaView to visualize the displacement field using the `viz/plot_dispwarp.py` Python script.
First, we start ParaView from the `examples/reverse-2d` directory.
Before running the `viz/plot_dispwarp.py` Python script as described in {ref}`sec-paraview-python-scripts`, we set the simulation name in the ParaView Python Shell.

```{code-block} python
---
caption: Set the simulation in the ParaView Python Shell.
---
>>> SIM = "step02_gravity_refstate"
```

:::{figure-md} fig:example:reverse:2d:step02:solution
<img src="figs/step02-solution.*" alt="Solution for Step 2. The colors indicate the magnitude of the displacement." width="100%"/>

Solution for Step 2.
The colors of the shaded surface indicate the magnitude of the displacement, which is zero.
The undeformed configuration is show by the gray wireframe.
:::
