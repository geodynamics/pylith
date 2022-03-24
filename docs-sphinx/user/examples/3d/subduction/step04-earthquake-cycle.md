# Step 4: Prescribed Earthquake Cycle

In Step 4, We combine the interseismic deformation in Step 3 with the coseismic slip in Step 2 to simulate two earthquake cycles.
We also include an earthquake on the splay fault. This illustrates how to include multiple earthquake sources on a single fault.
We use the same roller Dirichlet boundary conditions and combination of elastic and viscoelastic materials as we did in Step 3.

:::{figure-md} fig:example:subduction:3d:step04:diagram
<img src="figs/subduction3d_step04_diagram.*" alt="Diagram of Step 4 - A simple earthquake cycle combining the prescribed aseismic slip (creep) from Step 3 with prescribed coseismic slip for two earthquakes on the shallow portion of the subduction interface and one earthquake on the play fault. We impose roller Dirichlet boundary conditions on the lateral and bottom boundaries, except where they overlap with the slab and splay fault." width="100%"/>

Diagram of Step 4 - A simple earthquake cycle combining the prescribed aseismic slip (creep) from Step 3 with prescribed coseismic slip for two earthquakes on the shallow portion of the subduction interface and one earthquake on the play fault. We impose roller Dirichlet boundary conditions on the lateral and bottom boundaries, except where they overlap with the slab and splay fault.
:::

We create an array of three fault interfaces, one for the top of the slab (subduction interface), one for the bottom of the slab, and one for the splay fault.
The splay fault terminates into the fault on the top of the slab, so we must list the through-going fault on the top of the slab first.

```{code-block} cfg
---
caption: Excerpt from `step04.cfg`
---
# We prescribe slip on the top and bottom of the slab and on the splay fault.
[pylithapp.problem]
interfaces = [slab_top, slab_bottom, splay]

[pylithapp.problem.interfaces]
slab_top = pylith.faults.FaultCohesiveKin
slab_bottom = pylith.faults.FaultCohesiveKin
splay = pylith.faults.FaultCohesiveKin
```

:::{important}
When including intersecting faults, the through-going fault must be listed first in the array of fault interfaces. This ensures its cohesive cells are created before the adjacent fault that terminates into the through-going fault. For nonintersecting faults, the order in the list of fault interfaces does not matter.
:::

The settings for the fault interface on the bottom of the slab match those used in Step 3.
For the subduction interface, we want to impose creep on the deeper portion and earthquakes (coseismic slip) at specific times on the upper portion.
We create an array of earthquake sources, one for the creep and one for each of the earthquakes.
We want the earthquake to be imposed at specific times, so we set their origin time equal to the desire rupture time (100 years and 200 years) minus a value much smaller than the time step, so that roundoff errors do not result in the ruptures occurring one time step later than intended.
We use the same settings as we did in Step 3 for the creep earthquake source.
For the coseismic slip, we use a *SimpleGridDB* to impose a depth-dependent slip distribution that exactly complements the depth-dependent slip distribution of the creep.
Note that the slip time within an earthquake rupture is relative to the origin time, so we set the slip time to zero to coincide with the specified origin time.

```{code-block} cfg
---
caption: Excerpt from `step04.cfg`
---
# --- Skipping lines already discussed in Step 3 ---
eq_srcs = [creep, eq1, eq2]
eq_srcs.creep.origin_time = 0.0*year
eq_srcs.eq1.origin_time = 99.999*year ; 100*yr - small value
eq_srcs.eq2.origin_time = 199.999*year l 200*yr - small value

# Use the constant slip rate time function for the creep earthquake source.
eq_srcs.creep.slip_function = pylith.faults.ConstRateSlipFn

# Creep
[pylithapp.problem.interfaces.slab_top.eq_srcs.creep.slip_function]
slip_rate = spatialdata.spatialdb.SimpleGridDB
slip_rate.label = Slab top slip rate.
slip_rate.filename = spatialdb/fault_slabtop_creep.spatialdb
slip_rate.query_type = linear

slip_time = spatialdata.spatialdb.UniformDB
slip_time.label = Slip initiation time
slip_time.values = [slip-time]
slip_time.data = [0.0*year] ; Slip time is relative to origin time

# Earthquake 1
[pylithapp.problem.interfaces.slab_top.eq_srcs.eq1.slip_function]
slip = spatialdata.spatialdb.SimpleGridDB
slip.label = Slab top slip rate.
slip.filename = spatialdb/fault_slabtop_coseismic.spatialdb
slip.query_type = linear

slip_time = spatialdata.spatialdb.UniformDB
slip_time.label = Slip initiation time
slip_time.values = [slip-time]
slip_time.data = [0.0*year] ; Slip time is relative to origin time.

# Earthquake 2 (same as earthquake 1)
[pylithapp.problem.interfaces.slab_top.eq_srcs.eq2.slip_function]
slip = spatialdata.spatialdb.SimpleGridDB
slip.label = Slab top slip rate.
slip.filename = spatialdb/fault_slabtop_coseismic.spatialdb
slip.query_type = linear

slip_time = spatialdata.spatialdb.UniformDB
slip_time.label = Slip initiation time
slip_time.values = [slip-time]
slip_time.data = [0.0*year] ; Slip time is relative to origin time.
# --- Omitting output settings already discussed ---
```

The settings for the splay fault look very similar to those for the coseismic slip on the slab rupture patch in Step 2.
The primary difference is that we specify an origin time of 250 years.

```{code-block} cfg
---
caption: Excerpt from `step04.cfg`
---

[pylithapp.problem.interfaces.splay]
id = 102 ; id must be unique across all materials and faults
label = fault_splay ; Nodeset for the entire fault surface
edge = fault_splay_edge ; Nodeset for the buried edges

# We must define the quadrature information for fault cells.
# The fault cells are 2D (surface).
quadrature.cell = pylith.feassemble.FIATSimplex
quadrature.cell.dimension = 2

# Origin time for splay fault earthquake.
eq_srcs.rupture.origin_time = 249.999*year

# The slip time and final slip are defined in spatial databases.
[pylithapp.problem.interfaces.splay.eq_srcs.rupture.slip_function]
slip = spatialdata.spatialdb.UniformDB
slip.label = Splay fault slip.
slip.values = [left-lateral-slip, reverse-slip, fault-opening]
slip.data = [-1.0*m, 2.0*m, 0.0*m]

slip_time = spatialdata.spatialdb.UniformDB
slip_time.label = Slip initiation time
slip_time.values = [slip-time]
slip_time.data = [0.0*year] ; Relative to the origin time
# --- Omitting output settings already discussed ---
```

```{code-block} console
---
caption: Run Step 4 simulation
---
$ pylith step04.cfg mat_viscoelastic.cfg solver_fieldsplit.cfg
```

The simulation will produce sixteen pairs of HDF5/Xdmf files, beginning with `step04`, in the `output` directory:

**step04-domain.h5[.xmf]**  Time series of the solution field over the domain.

**step04-groundsurf.h5[.xmf]**  Time series of the solution field over the ground surface.

**step04-slab_info.h5[.xmf]**  Properties for the slab material.

**step04-slab.h5[.xmf]**  Time series of the state variables (stress and strain) for the slab material.

**step04-wedge_info.h5[.xmf]**  Properties for the wedge material.

**step04-wedge.h5[.xmf]**  Time series of the state variables (stress and strain) for the wedge material.

**step04-crust_info.h5[.xmf]**  Properties for the crust material.

**step04-crust.h5[.xmf]**  Time series of the tate variables (stress and strain) for the crust material.

**step04-mantle_info.h5[.xmf]**  Properties for the mantle material.

**step04-mantle.h5[.xmf]**  Time series of the state variables (stress and strain) for the mantle material.

**step04-fault-slabbot_info.h5[.xmf]**  Fault orientation and rupture information for the bottom of the slab.

**step04-fault-slabbot.h5[.xmf]**  Time series of slip and traction changes for the bottom of the slab.

**step04-fault-slabtop_info.h5[.xmf]**  Fault orientation and rupture information for the top of the slab.

**step04-fault-slabtop.h5[.xmf]**  Time series of slip and traction changes for the top of the slab.

**step04-fault-splay_info.h5[.xmf]**  Fault orientation and rupture information for the splay fault.

**step04-fault-splay.h5[.xmf]**  Time series of slip and traction changes for the splay fault.

{numref}`fig:example:subduction:3d:step04`, which was created using the ParaView Python script `plot_dispwarp.py`, shows the deformation exaggerated by a factor of 5,000 at the final time step of t=300*yr.
Compared to the solution in Step 3, we see the earthquakes have reduced the deformation in the crust and accretionary wedge.

:::{figure-md} fig:example:subduction:3d:step04
<img src="figs/subduction3d_step04_soln.*" alt="Solution over the domain for Step 4 at {math}`t=300`yr. The colors indicate the z-displacement and we have exaggerated the deformation by a factor of 5,000." width="100%"/>

Solution over the domain for Step 4 at {math}`t=300`yr. The colors indicate the z-displacement and we have exaggerated the deformation by a factor of 5,000.
:::

## Exercises

* Adjust the timing of the earthquake rupture sequence. How does this affect the deformation?
* Add additional earthquakes with different depth variations in slip, keeping the total equal to the overall slip rate.
* Adjust the nodesets in CUBIT/Trelis so that the splay fault and the deeper portion of the subduction interface form the through-going fault and the upper portion of the subduction interface is the secondary fault. How does this affect the stress accumulation in the crust and upper mantle?
