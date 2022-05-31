# Step 7: Slip on Two Faults and Maxwell Viscoelastic Materials

In this example we replace the linear, isotropic elastic bulk rheology in the slab with a linear, isotropic Maxwell viscoelastic rheology.
We also switch from a static simulation to a quasistatic simulation to compute the time-dependent relaxation in the slab.
We use the same boundary conditions as in Step 6.

% Metadata extracted from parameter files.
```{include} step07_twofaults_maxwell-synopsis.md
```

## Simulation parameters

The parameters specific to this example are in `step07_twofaults_maxwell.cfg` and include:

* `pylithapp.metadata` Metadata for this simulation. Even when the author and version are the same for all simulations in a directory, we prefer to keep that metadata in each simulation file as a reminder to keep it up-to-date for each simulation.
* `pylithapp` Parameters defining where to write the output.
* `pylithapp.problem` Parameters for the time stepping and solution field with displacement and Lagrange multiplier subfields.
* `pylithapp.problem.material` Parameters for the linear Maxwell viscoelastic bulk rheology for the slab.
* `pylithapp.problem.fault` Parameters for prescribed slip on the two faults.

We use a very short relaxation time of 20 years, so we run the simulation for 100 years with a time step of 4 years.
We use a starting time of -4 years so that the first time step will advance the solution time to 0 years.

:::{important}
Both faults contain one end that is buried within the domain.
The splay fault ends where it meets the main fault.
When PyLith inserts cohesive cells into a mesh with buried edges (in this case a point), we must identify these buried edges so that PyLith properly adjusts the topology along these edges.

For properly topology of the cohesive cells, the main fault _must_ be listed first in the array of faults so that it will be created before the splay fault.
:::

```{code-block} console
---
caption: Run Step 7 simulation
---
$ pylith step07_twofaults_maxwell.cfg

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

25 TS dt 0.2 time 4.8
    0 SNES Function norm 2.589196279152e-05 
    Linear solve converged due to CONVERGED_ATOL iterations 339
    1 SNES Function norm 6.643244443329e-13 
  Nonlinear solve converged due to CONVERGED_FNORM_ABS iterations 1
26 TS dt 0.2 time 5.
 >> /Users/baagaard/software/unix/py39-venv/pylith-debug/lib/python3.9/site-packages/pylith/problems/Problem.py:201:finalize
 -- timedependent(info)
 -- Finalizing problem.
```

From the end of the output written to the terminal window, we see that the simulation advanced the solution 26 time steps.
The PETSc TS display time in the nondimensional units, so a time of 5 corresponds to 100 years.

## Visualizing the results

The `output` directory contains the simulation output.
Each "observer" writes its own set of files, so the solution over the domain is in one set of files, the boundary condition information is in another set of files, and the material information is in yet another set of files.
The HDF5 (`.h5`) files contain the mesh geometry and topology information along with the solution fields.
The Xdmf (`.xmf`) files contain metadata that allow visualization tools like ParaView to know where to find the information in the HDF5 files.
To visualize the data using ParaView or Visit, load the Xdmf files.

In {numref}`fig:example:reverse:2d:step07:solution` we use ParaView to visualize the y displacement field using the `viz/plot_dispwarp.py` Python script.
First, we start ParaView from the `examples/reverse-2d` directory.
Before running the `viz/plot_dispwarp.py` Python script as described in {ref}`sec-paraview-python-scripts`, we set the simulation name in the ParaView Python Shell.

```{code-block} python
---
caption: Set the simulation and exaggeration in the ParaView Python Shell.
---
>>> SIM = "step07_twofaults_maxwell"
>>> FIELD_COMPONENT = "X"
```

:::{figure-md} fig:example:reverse:2d:step07:solution
<img src="figs/step07-solution.*" alt="Solution for Step 7. The colors indicate the magnitude of the x displacement, and the deformation is exaggerated by a factor of 1000." width="100%"/>

Solution for Step 7 at t=100 years.
The colors of the shaded surface indicate the magnitude of the x displacement, and the deformation is exaggerated by a factor of 1000.
The undeformed configuration is show by the gray wireframe.
Viscoelastic relaxation results in significant deformation in the `slab` material.
:::
