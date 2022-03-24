(sec:example:subduction:3d:step03)=
# Step 3: Prescribed Aseismic Creep and Interseismic Deformation

We now increase the complexity of our fault model by simulating the interseismic deformation associated with the subducting slab.
We approximate the motion of the Juan de Fuca Plate subducting under the North American Plate by introducing aseismic slip (creep) on the bottom of the slab and the deeper portion of the subduction interface; we keep the interface between the subduction interface and the accretionary wedge and shallow crust locked.
As in Step 2, we will use the linear Maxwell viscoelastic constitutive model for the slab and mantle.
{numref}`fig:example:subduction:3d:step03:diagram` summarizes the problem description.

:::{figure-md} fig:example:subduction:3d:step03:diagram
<img src="figs/subduction3d_step03_diagram.*" alt="Diagram of Step 3 - Prescribed aseismic slip (creep) and interseismic deformation for the subducting slab. We prescribe steady, uniform creep on the bottom of the slab and deeper portion of the subduction interface. We impose roller Dirichlet boundary conditions on the lateral and bottom boundaries, except where they overlap with the slab and splay fault." width="100%"/>

Diagram of Step 3 - Prescribed aseismic slip (creep) and interseismic deformation for the subducting slab. We prescribe steady, uniform creep on the bottom of the slab and deeper portion of the subduction interface. We impose roller Dirichlet boundary conditions on the lateral and bottom boundaries, except where they overlap with the slab and splay fault.
:::

% Fault
With slip on the top and bottom of the slab, our fault interfaces array contains two components, one for the top of the slab (subduction interface), **slab_top**, and one for the bottom of the slab, **slab_bottom**.
We use the *FaultCohesiveKin* object for each of these interfaces since we want to prescribe the slip.

```{code-block} cfg
---
caption: Excerpt from `step03.cfg`
---
[pylithapp.problem]
interfaces = [slab_bottom, slab_top]

[pylithapp.problem.interfaces]
slab_bottom = pylith.faults.FaultCohesiveKin
slab_top = pylith.faults.FaultCohesiveKin
```

We specify the **id** used to identify the cohesive cells for this fault so that it is unique among all materials and faults.
We also specify the appropriate nodesets identifying the entire fault surface and the buried edges.
Some portions of the bottom of the slab are perfectly horizontal, so our procedure that uses the vertical direction and the fault normal to set the along-strike and up-dip shear components breaks down.
We remedy this by tweaking the **up_dir** direction from being completely vertical (0,0,1) to tilting slightly to the west.
This results in consistent along-strike and up-dip directions across the fault surface.
For the aseismic slip we use a constant slip rate time function (*ConstRateSlipFn*) with *UniformDB* spatial databases to specify the constant, uniform oblique slip rate of 2.0 cm/yr of left-lateral motion and 4.0 cm/yr of normal motion.
Note that slip on the bottom of the subducting slab has the opposite sense of motion as that on the top of the slab.

```{code-block} cfg
---
caption: Excerpt from `step03.cfg`
---
[pylithapp.problem.interfaces.slab_bottom]
id = 100 ; Must be different from ids used for materials
label = fault_slabbot ; Nodeset for the entire fault surface
edge = fault_slabbot_edge ; Nodeset for the buried edges
# Give slight westward tilt to the up_dir to avoid ambiguous
# directions for the shear components on the horizontal portions of the
# fault.
up_dir = [-0.1,0,0.9]

# We must define the quadrature information for fault cells.
# The fault cells are 2D (surface).
quadrature.cell = pylith.feassemble.FIATSimplex
quadrature.cell.dimension = 2

# Use the constant slip rate time function.
eq_srcs.rupture.slip_function = pylith.faults.ConstRateSlipFn

# The slip time and final slip are defined in spatial databases.
[pylithapp.problem.interfaces.slab_bottom.eq_srcs.rupture.slip_function]
slip_rate = spatialdata.spatialdb.UniformDB
slip_rate.label = Slab bottom slip rate.
slip_rate.values = [left-lateral-slip, reverse-slip, fault-opening]
slip_rate.data = [+2.0*cm/year, -4.0*cm/year, 0.0*cm/year]

slip_time = spatialdata.spatialdb.UniformDB
slip_time.label  = Slip initiation time
slip_time.values = [slip-time]
slip_time.data = [0.0*year]

[pylithapp.problem.interfaces.slab_bottom.output]
writer = pylith.meshio.DataWriterHDF5
writer.filename = output/step03-fault-slabbot.h5
vertex_info_fields = [normal_dir, strike_dir, dip_dir]
```

The parameters for the top of the slab (subduction interface) closely resemble those for the bottom of the slab.
The main difference is that we use a *SimpleGridDB* to define a depth variation in the slip rate.
The fault is locked at depths above 45 km and increases linearly to the same slip rate as the bottom of the slab at a depth of 60 km.

```{code-block} cfg
---
caption: Excerpt from `step03.cfg`
---
[pylithapp.problem.interfaces.slab_top]
id = 101 ; Must be different from ids used for materials
label = fault_slabtop ; Nodeset for the entire fault surface
edge = fault_slabtop_edge ; Nodeset for the buried edges

# We must define the quadrature information for fault cells.
# The fault cells are 2D (surface).
quadrature.cell = pylith.feassemble.FIATSimplex
quadrature.cell.dimension = 2

# Use the constant slip rate time function.
eq_srcs.rupture.slip_function = pylith.faults.ConstRateSlipFn

# The slip time and final slip are defined in spatial databases.
[pylithapp.problem.interfaces.slab_top.eq_srcs.rupture.slip_function]
slip_rate = spatialdata.spatialdb.SimpleGridDB
slip_rate.label = Slab top slip rate.
slip_rate.filename = spatialdb/fault_slabtop_creep.spatialdb
slip_rate.query_type = linear

slip_time = spatialdata.spatialdb.UniformDB
slip_time.label  = Slip initiation time
slip_time.values = [slip-time]
slip_time.data = [0.0*year]

[pylithapp.problem.interfaces.slab_top.output]
writer = pylith.meshio.DataWriterHDF5
writer.filename = output/step03-fault-slabtop.h5
vertex_info_fields = [normal_dir, strike_dir, dip_dir]
```

We do not want the boundaries to constrain the motion of the subducting slab, so we use the nodesets that exclude vertices on the subducting slab.
Furthermore, PyLith does not permit overlap between the fault interfaces and Dirichlet boundary conditions.
This is why we exclude vertices on the splay fault in these nodesets as well.
We only update the name of the nodeset for the -x, -y, and +y boundaries.

```{code-block} cfg
---
caption: Excerpt from `step03.cfg`
---
# -x face
[pylithapp.problem.bc.x_neg]
label = boundary_xneg_noslab

# -y face
[pylithapp.problem.bc.y_neg]
label = boundary_yneg_noslab

# +y face
[pylithapp.problem.bc.y_pos]
label = boundary_ypos_noslab
```

```{code-block} console
---
caption: Run Step 3 simulation
---
$ pylith step03.cfg mat_viscoelastic.cfg solver_fieldsplit.cfg

```

The simulation will produce fourteen pairs of HDF5/Xdmf files, beginning with `step03`, in the `output` directory:

**step03-domain.h5[.xmf]**  Time series of the solution field over the domain.

**step03-groundsurf.h5[.xmf]**  Time series of the solution field over the ground surface.

**step03-slab_info.h5[.xmf]**  Properties for the slab material.

**step03-slab.h5[.xmf]**  Time series of the state variables (stress and strain) for the slab material.

**step03-wedge_info.h5[.xmf]**  Properties for the wedge material.

**step03-wedge.h5[.xmf]**  Time series of the state variables (stress and strain) for the wedge material.

**step03-crust_info.h5[.xmf]**  Properties for the crust material.

**step03-crust.h5[.xmf]**  Time series of the tate variables (stress and strain) for the crust material.

**step03-mantle_info.h5[.xmf]**  Properties for the mantle material.

**step03-mantle.h5[.xmf]**  Time series of the state variables (stress and strain) for the mantle
material.

**step03-fault-slabbot_info.h5[.xmf]**  Fault orientation and rupture information for the bottom of the slab.

**step03-fault-slabbot.h5[.xmf]**  Time series of slip and traction changes for the bottom of the slab.

**step03-fault-slabtop_info.h5[.xmf]**  Fault orientation and rupture information for the top of the slab.

**step03-fault-slabtop.h5[.xmf]**  Time series of slip and traction changes for the top of the slab.

As in Step 2, there are two pairs of HDF5/Xdmf files for each fault; one set for the fault orientation and rupture information and one set for the time series of slip and change in tractions.

{numref}`fig:example:subduction:3d:step03`, which was created using the ParaView Python script `plot_dispwarp.py`, shows the deformation exaggerated by a factor of 5,000 at the final time step of t=200*yr.
Notice that there are some local edge effects associated with the unconstrained degrees of freedom at the intersection of the boundaries and fault surfaces.

:::{figure-md} fig:example:subduction:3d:step03
<img src="figs/subduction3d_step03_soln.*" alt="Solution over the domain for Step 2 at {math}`t=200`yr. The colors indicate the x-displacement and we have exaggerated the deformation by a factor of 5,000." width="100%"/>

Solution over the domain for Step 2 at {math}`t=200`yr. The colors indicate the x-displacement and we have exaggerated the deformation by a factor of 5,000.
:::

## Exercises

* Adjust the locking depth for the subduction interface. How does this affect the spatial distribution of the change in tractions on the fault interfaces?
* Increase the rigidity of the slab and decrease the rigidity of the wedge and/or crust. How do these affect the change in tractions on the fault interfaces?
