# Examples: 3-D Subduction Zone

This suite of examples demonstrates use of a wide variety of features
and the general workflow often used in research simulations.  The
examples consistent of a step-by-step sequence of eight problems
involving a 3-D subduction zone. They focus on modeling the
deformation associated with the the subducting slab, including
interseismic deformation with aseismic slip (creep) and viscoelastic
relaxation, coseismic slip on the slab interface and a splay fault,
and slow slip events on the subduction interface. We want to account
for the 3-D material properties associated with different elastic
properties for the subducting slab, mantle, continental crust, and an
accretionary wedge. All of the examples use the same mesh, generated
using CUBIT/Trelis.

Directory structure:
* **mesh**: mesh related files
* **output**: simulation output files, created automatically by PyLith
* **spatialdb**: spatial and temporal database related files
* **viz**: ParaView Python scripts

Examples:

* `Step 1`: Axial compression
* `Step 2`: Prescribed coseismic slip and postseismic relaxation
* `Step 3`: Prescribed aseismic creep and interseismic deformation
* `Step 4`: Prescribed earthquake cycle
* `Step 5`: Spontaneous rupture driven by subducting slab (not yet working for v3.0)
* `Step 6`: Prescribed slow-slip event
* `Step 7a,b`: Inversion of slow-slip event using 3-D Green's functions (not yet working for v3.0)
* `Step 8a,b,c`: Stress field due to gravitational body forces


## Mesh generation using CUBIT/Trelis (optional)

We use CUBIT/Trelis to generate the finite-element mesh. Due to its
size, we do not include the finite-element mesh in the PyLith source
or binary distributions. If you do not have CUBIT/Trelis, you can
download the mesh from
https://wiki.geodynamics.org/software:pylith:examples:files and skip
generating the mesh.

```
# Generate generate_surfs.jou
# Make sure you are in the 'mesh' directory and then run the Python
# script to generate the journal file 'generate_surfs.jou'.
./generate_surfjou.py
```

The next step is to use CUBIT/Trelis to run the `generate_surfs.jou`
journal file to generate the spline surfaces for the slab and splay
fault and save them as ACIS surfaces.

After you generate the ACIS surface files, run the `mesh_tet.jou`
journal file to construct the geometry, and generate the mesh. In the
end you will have an Exodus-II file `mesh_tet.exo`, which is a NetCDF
file, in the `mesh` directory.


## Step01: Axial Compression

We start with a very simple example of axial compression in the
east-west direction with purely elastic material properties, and no
faults. 

To run this example:
```
pylith step01_axialdisp.cfg mat_elastic.cfg
```


## Step02: Prescribed Coseismic Slip and Postseismic Relaxation

In this example we model the postseismic relaxation of the deep slab
and mantle resulting from coseismic slip on a fault patch in the
central portion of the subduction (top of the slab) interface. 

To run this example:

```
pylith step02_coseismic.cfg mat_viscoelastic.cfg solver_fieldsplit.cfg
```


## Step03: Prescribed Aseismic Creep and Interseismic Deformation

We now increase the complexity of our fault model by simulating the
interseismic deformation associated with the subducting slab.

To run this example:
```
pylith step03_interseismic.cfg mat_viscoelastic.cfg solver_fieldsplit.cfg
```


## Step04: Prescribed Earthquake Cycle

In Step 4 we combine the interseismic deformation in Step 3 with the
coseismic slip in Step 2 to simulate two earthquake cycles. We also
include an earthquake on the splay fault.

To run this example:
```
pylith step04_eqcycle.cfg mat_viscoelastic.cfg solver_fieldsplit.cfg
```

## Step05: Spontaneous Rupture Driving by Subducting Slab

This example is not yet complete.

## Step06: Prescribed Slow-Slip Event

This example simulates a simple slow slip event (SSE) on the
subduction interface, in which the entire patch slips simultaneously
with an amplitude that grows with time. 

To run this example:
```
# Generate the spatial database files: 
#    - fault_slabtop_slowslip.spatialdb
#    - fault_slabtop_slowslip.timedb
cd spatialdb && ./generate_slowslip.py
# Change back to the subduction directory and run PyLith
$ cd ..
$ pylith step06_slowslip.cfg mat_elastic.cfg solver_fieldsplit.cfg
```


## Step07: Inversion of Slow-Slip Event using 3-D Green's Functions

**WARNING**: This example will not work with PyLith v3.0. It uses
Green's functions, which has not been implemented in v3.0. Use
the PyLith v2.2.1 release to generate static Green's functions.

This example is a three-dimensional analog of generating Green's
functions in 2-D and is a more realistic example of how PyLith can be
used to perform geodetic inversions. We divide generating Green's
functions for slip impulses on the central rupture patch of the
subduction interface two sub-problems:

  * **Step 7a**: Left-lateral slip component
  * **Step 7b**: Reverse slip component

Although PyLith can generate the two components in one simulation, we
often prefer to speed up the process by running simulations for each
of the components at the same time using multiple processes on a
cluster.

To run this example:
```
pylith --problem=pylith.problems.GreensFns step07a_leftlateral.cfg mat_elastic.cfg solver_fieldsplit.cfg
pylith --problem=pylith.problems.GreensFns step07b_reverse.cfg mat_elastic.cfg solver_fieldsplit.cfg
```

Before we can run the inversion, we post-process the output from Step
6 to create synthetic data. 
```
# Generate synthetic GPS data
./make_synthetic_gpsdisp.py
```

We perform a simple inversion using the `slip_invert.py` script,
with parameters defined in `slip_invert.cfg`.
```
./slip_invert.py
```


## Step08: Stress Field Due to Gravitational Body Forces

This example demonstrates the use of gravitational body forces as well
as the use of initial stresses to balance the body forces. This
involves enabling gravity within our domain with Dirichlet roller
boundary conditions on the lateral and bottom boundaries; we do not
include faults in this example.  We also demonstrate what happens when
the initial stresses are not in balance with the gravitational
stresses, and show how viscoelastic problems with gravitational
stresses will in general not reach a steady-state solution. The
example is divided into three sub-problems:

  * **Step 8a**: Gravitational body forces with 3-D density variations
    in elastic materials and initial stresses for a uniform density.
  * **Step 8b**: Gravitational body forces with 3-D density variations
    in incompressible elastic materials.
  * **Step 8c**: Gravitational body forces with 3-D density variations
    in elastic and viscoelastic materials and initial stresses from
    Step 8a plus finite strain formulation (does not reach a
    steady-state solution).

To run this example:
```
# Step 8a
pylith step08a_gravity_refstate.cfg mat_elastic.cfg solver_algebraicmultigrid.cfg

# Step 8b
pylith step08b_gravity_incompressible.cfg mat_incompressible.cfg solver_algebraicmultigrid.cfg

# Step 8c
pylith step08c_gravity_viscoelastic.cfg mat_viscoelastic.cfg solver_algebraicmultigrid.cfg
```

\end{shell}

