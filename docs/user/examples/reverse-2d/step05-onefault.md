# Step 5: Static Coseismic Slip

% Metadata extracted from parameter files.
```{include} step05_onefault-synopsis.md
```

## Simulation parameters

This example involves a static simulation that solves for the deformation from prescribed coseismic slip on the main fault.
We specify 2 meters of reverse slip.
{numref}`fig:example:reverse:2d:step05:diagram` shows the boundary conditions on the domain.
The parameters specific to this example are in `step05_onefault.cfg`.

:::{figure-md} fig:example:reverse:2d:step05:diagram
<img src="figs/step05-diagram.*" alt="" scale="75%">

Boundary conditions for static coseismic slip on the main fault.
We prescribe 2 meters of reverse slip with roller boundary conditions on the lateral sides and bottom of the domain.
:::

:::{important}
In 2D simulations slip is specified in terms of opening and left-lateral components.
This provides a consistent, unique sense of slip that is independent of the fault orientation.
For our geometry in this example, right lateral slip corresponds to reverse slip on the dipping fault.
:::

:::{important}
The main fault contains one end that is buried within the domain.
When PyLith inserts cohesive cells into a mesh with buried edges (in this case a point), we must identify these buried edges so that PyLith properly adjusts the topology along these edges.
:::

We adjust the solution field to include both displacement and the Lagrange multiplier associated with the fault.
For uniform prescribed slip we use a `UniformDB`.

```{code-block} cfg
---
caption: Parameters for earthquake rupture on the main reverse fault in Step 5.
---
[pylithapp.problem]
solution = pylith.problems.SolnDispLagrange

[pylithapp.problem.interfaces.fault]
label = fault
label_value = 20
edge = fault_end
edge_value = 21
observers.observer.data_fields = [slip]

[pylithapp.problem.interfaces.fault.eq_ruptures.rupture]
db_auxiliary_field = spatialdata.spatialdb.UniformDB
db_auxiliary_field.description = Fault rupture auxiliary field spatial database
db_auxiliary_field.values = [initiation_time, final_slip_left_lateral, final_slip_opening]
db_auxiliary_field.data = [0.0*s, -2.0*m, 0.0*m]
```

## Running the simulation

```{code-block} console
---
caption: Run Step 5 simulation
---
$ pylith step05_onefault.cfg

# The output should look something like the following.
 >> /software/unix/py3.12-venv/pylith-debug/lib/python3.12/site-packages/pylith/apps/PyLithApp.py:77:main
 -- pylithapp(info)
 -- Running on 1 process(es).
 >> /software/unix/py3.12-venv/pylith-debug/lib/python3.12/site-packages/pylith/meshio/MeshIOObj.py:38:read
 -- meshiopetsc(info)
 -- Reading finite-element mesh
 >> /src/cig/pylith/libsrc/pylith/meshio/MeshIO.cc:85:void pylith::meshio::MeshIO::read(pylith::topology::Mesh *, const bool)
 -- meshiopetsc(info)
 -- Component 'reader': Domain bounding box:
    (-100000, 100000)
    (-100000, 0)

# -- many lines omitted --

 >> /software/unix/py3.12-venv/pylith-debug/lib/python3.12/site-packages/pylith/problems/TimeDependent.py:132:run
 -- timedependent(info)
 -- Solving problem.
0 TS dt 0.01 time 0.
    0 SNES Function norm 2.129295960330e-02
    Linear solve converged due to CONVERGED_ATOL iterations 30
    1 SNES Function norm 2.046522086472e-12
  Nonlinear solve converged due to CONVERGED_FNORM_ABS iterations 1
1 TS dt 0.01 time 0.01
 >> /software/unix/py3.12-venv/pylith-debug/lib/python3.12/site-packages/pylith/problems/Problem.py:199:finalize
 -- timedependent(info)
 -- Finalizing problem.
```

From the end of the output written to the terminal window, we see that the linear solver converged in 30 iterations and met the absolute convergence tolerance (`ksp_atol`).
As we expect for this linear problem, the nonlinear solver converged in 1 iteration.

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
