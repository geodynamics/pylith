# Step 6: Prescribed Slow-Slip Event

This example simulates a simple slow slip event (SSE) on the subduction interface, in which the entire patch slips simultaneously with an amplitude that grows with time.
We impose a constant rake angle of 110 degrees, and a time duration of 30 days.
The time duration is much shorter than the Maxwell time for our viscoelastic materials, so we use elastic material properties (as we did in Step 1).

:::{figure-md} fig:example:subduction:3d:step06:diagram
<img src="figs/subduction3d_step06_diagram.*" alt="Diagram of Step 6 - Prescribed slow-slip event on the subduction interface. This quasistatic simulation prescribes a Gaussian slip distribution on the central rupture patch of the subduction interface, purely elastic material properties, and roller boundary conditions on the lateral (north, south, east, and west) and bottom boundaries." width="100%"/>

Diagram of Step 6 - Prescribed slow-slip event on the subduction interface. This quasistatic simulation prescribes a Gaussian slip distribution on the central rupture patch of the subduction interface, purely elastic material properties, and roller boundary conditions on the lateral (north, south, east, and west) and bottom boundaries.
:::

The only time dependence in this problem is the time evolution of slip, so we set the duration of the simulation to match the duration of the slow slip event.
We use a time step of 2.0 days to insure that we resolve the temporal evolution of the slip.

```{code-block} cfg
---
caption: Excerpt from `step06.cfg`
---
[pylithapp.problem.formulation.time_step]
total_time = 30.0*day
dt = 2.0*day
```

The results in this example will be used to simulate output at fake continuous GPS (cGPS) stations in Step 7, so we add an output manager for saving the solution at specific points (`OutputSolnPoints`) in addition to our output managers over the domain and top surface:

```{code-block} cfg
---
caption: Excerpt from `step06.cfg`
---
[pylithapp.problem.implicit]
output = [domain, subdomain, cgps_sites]

# Default output is for the entire domain.
# We need to set the type of output for the subdomain and points.
output.subdomain = pylith.meshio.OutputSolnSubset
output.cgps_sites = pylith.meshio.OutputSolnPoints
```

For the point output we specify the output data writer, the file containing the list of cGPS stations and the coordinate system associated with the station locations.
The format of the station file is whitespace separated columns of station name and then the coordinates of the station.
See {numref}`sec:format:PointsList` for more information.

```{code-block} cfg
---
caption: Excerpt from `step06.cfg`
---
[pylithapp.problem.formulation.output.cgps_sites]
writer = pylith.meshio.DataWriterHDF5
writer.filename = output/step06-cgps_sites.h5

# File with coordinates of cGPS stations.
reader.filename = cgps_sites.txt

# Specify coordinate system used in cGPS station file.
coordsys = spatialdata.geocoords.CSGeo
coordsys.space_dim = 3
coordsys.datum_horiz = WGS84
coordsys.datum_vert = mean sea level
```

The fault parameters are very similar to those in Step 2, in which we also prescribed slip on the subduction interface patch.
The primary difference is that we use a user-defined slip time history function (*TimeHistorySlipFn*).
This slip time function requires spatial databases for the amplitude of the final slip and slip initiation time, and a time history file specifying the normalized amplitude as a function of time.
Additionally, to illustrate PyLith's ability to use spatial databases with points in other, but compatible, georeferenced coordinate systems, we specify the slip distribution using geographic (longitude and latitude) coordinates.

```{code-block} cfg
---
caption: Excerpt from `step06.cfg`
---
[pylithapp.problem]
# We prescribe slip on the slab fault patch.
interfaces = [slab]

[pylithapp.problem.interfaces]
slab = pylith.faults.FaultCohesiveKin

[pylithapp.problem.interfaces.slab]
# Nodeset corresponding to the fault patch and buried edge.
label = fault_slabtop_patch
edge = fault_slabtop_patch_edge

# We must define the quadrature information for fault cells.
# The fault cells are 2D (surface).
quadrature.cell = pylith.feassemble.FIATSimplex
quadrature.cell.dimension = 2

# We use a time history slip function.
[pylithapp.problem.interfaces.slab.eq_srcs.rupture]
slip_function = pylith.faults.TimeHistorySlipFn

# The slip is defined in a SimpleGridDB spatial database with linear interpolation.
[pylithapp.problem.interfaces.slab.eq_srcs.rupture.slip_function]
slip = spatialdata.spatialdb.SimpleGridDB
slip.label = Gaussian slip distribution for SSE
slip.filename = spatialdb/fault_slabtop_slowslip.spatialdb
slip.query_type = linear

# We use a UniformDB to specify the slip initiation time.
slip_time = spatialdata.spatialdb.UniformDB
slip_time.label = Slip initiation time
slip_time.values = [slip-time]
slip_time.data = [0.0*year]

# We use a temporal database to provide the slip time history.
time_history.label = Time history of slip
time_history.filename = spatialdb/fault_slabtop_slowslip.timedb
```

You will notice that the `spatialdb` directory does not contain the `fault_slabtop_slowslip.spatialdb` and `fault_slabtop_slowslip.timedb` files.
We use the `generate_slowslip.py` Python script to generate these files as an illustration of how to use Python to generate more simple spatial variations and the *SimpleGridAscii* object to write spatial database files.
This script reads parameters from `generate_slowslip.cfg` to generate a Gaussian slip distribution in geographic coordinates, along with a temporal database providing the slip amplitudes at different times.

:::{tip}
The `generate_slowslip.py` script is one of several examples where we make use of the Python interface to the spatialdata package.
This provides useful methods for handling coordinate systems and spatial databases.
:::

To run the simulation, first run the Python script to generate the spatial database files, and then run PyLith.

```{code-block} console
---
caption: Run Step 6 simulation
---
# Generate the spatial database files
$ cd spatialdb && ./generate_slowslip.py
$ ls fault_slabtop_slowslip.*
# You should see
fault_slabtop_slowslip.spatialdb  fault_slabtop_slowslip.timedb
# Change back to the subduction directory and run PyLith
$ cd ..
$ pylith step06.cfg mat_elastic.cfg solver_fieldsplit.cfg
```

The problem will produce thirteen pairs of HDF5/Xdmf files:

**step06-domain.h5[.xmf]**  Time series of the solution field over the domain.

**step06-groundsurf.h5[.xmf]**  Time series of the solution field over the ground surface.

**step06-cgps_sites.h5[.xmf]**  Time series of the solution field at the cGPS sites.

**step06-slab_info.h5[.xmf]**  Properties for the slab material.

**step06-slab.h5[.xmf]**  Time series of the state variables (stress and strain) for the slab material.

**step06-wedge_info.h5[.xmf]**  Properties for the wedge material.

**step06-wedge.h5[.xmf]**  Time series of the state variables (stress and strain) for the wedge material.

**step06-crust_info.h5[.xmf]**  Properties for the crust material.

**step06-crust.h5[.xmf]**  Time series of the tate variables (stress and strain) for the crust material.

**step06-mantle_info.h5[.xmf]**  Properties for the mantle material.

**step06-mantle.h5[.xmf]**  Time series of the state variables (stress and strain) for the mantle material.

**step06-fault-slab_info.h5[.xmf]**  Fault orientation and rupture information for the top of the slab.

**step06-fault-slab.h5[.xmf]**  Time series of slip and traction changes for the top of the slab.

The additional HDF5 file that was not present in previous examples is `step06-cgps_sites.h5`, which contains the displacements at the fake cGPS sites.

{numref}`fig:example:subduction:3d:step06`, which was created using ParaView, shows the surface vertical displacement along with horizontal displacement vectors at the cGPS sites, superimposed on contours of the applied slip at {math}`t = 24` days.

:::{figure-md} fig:example:subduction:3d:step06
<img src="figs/subduction3d_step06_soln.*" alt="Solution for Step 6. The colors indicate the vertical displacement, the vectors represent the horizontal displacements at fake cGPS sites, and the contours represent the applied slip at {math}`t = 24` days." width="100%"/>

Solution for Step 6. The colors indicate the vertical displacement, the vectors represent the horizontal displacements at fake cGPS sites, and the contours represent the applied slip at {math}`t = 24` days.
:::

## Exercises

* Change spatial distribution and time history of slip.
  * Edit `generate_slowslip.cfg` to change spatial and temporal distributions, and edit `step06.cfg` to change the time duration and/or time step size.
* Add propagation of the slow slip (spatial variation of slip initiation time).
  * Either alter Python script to produce a spatial database of slip initiation times, or write a new script. Can you produce a more realistic-looking slow slip event?
