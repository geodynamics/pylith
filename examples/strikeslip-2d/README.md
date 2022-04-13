# Examples: 2-D Strike-Slip Fault

This suite of examples demonstrates some basic concepts of using
PyLith with a single throughgoing strike-slip fault and a very simple
meshing example.
Concepts common to all of the steps include:

* Dirichlet boundary conditions
* Isotropic, linear elasticity with two materials and a contrast in shear
  modulus across the fault
* UniformDB spatial database for specifying values for material properties
* Kinematic fault slip
* Output of the solution over the domain, boundaries, fault, and materials
* Output of auxiliary information for boundary conditions and
  materials

## Meshing: Meshing a 2-D geometry including a throughgoing fault

We create the mesh using CUBIT/Trelis; however the resulting mesh is
provided, so you can skip creating the mesh if you do not have
CUBIT/Trelis. The associated journal files are: `geometry.jou`,
`createbc.jou`, `gradient.jou`, `mesh_tri.jou` (triangule cells) and
`mesh_quad.jou` (quadrilateral mesh).

## Step01: Prescribed slip with fixed Dirichlet boundary conditions

Zero-displacement Dirichlet boundary conditions on the +x, and -x
boundaries, and a single coseismic slip event. Features used in this
simulation include:

* Static simulation
* ZeroDB spatial database for specifying values for boundary conditions
* UniformDB spatial database for specifying values for fault slip and
  slip initiation time
* Use of precsribed static fault slip

The simulation parameters are in the `pylithapp.cfg` and
`step01_slip.cfg` files.

To run the example:
```
pylith step01_slip.cfg
```

## Step02: Prescribed slip with time-dependent Dirichlet boundary conditions

Velocity boundary conditions applied in the y-direction, and a single
earthquake rupture with slip matching the accumulated slip
deficit. Features used in this simulation include:

* Quasi-static, time-dependent simulation
* SimpleDB spatial database for specifying values for boundary conditions
* UniformDB spatial database for specifying values for fault slip and
  slip initiation time
* Rate conditions used to provide velocities on left and right sides

The simulation parameters are in the `pylithapp.cfg` and
`step02_slip_velbc.cfg` files.

To run the example:
```
pylith step02_slip_velbc.cfg
```

## Step03: Multiple ruptures with time-dependent Dirichlet boundary conditions

The same as problem Step02, but in this case we run the problem for twice as
long and have ruptures at 100 and 200 years, with the slip amount matching the
accumulated slip deficit at both times. Features used in this simulation
include:

* Quasi-static, time-dependent simulation
* SimpleDB spatial database for specifying values for boundary conditions
* UniformDB spatial database for specifying values for fault slip and
  slip initiation time
* Rate conditions used to provide velocities on left and right sides
* Use of multiple rupture sources

The simulation parameters are in the `pylithapp.cfg` and
`step03_multislip_velbc.cfg` files.

To run the example:
```
pylith step03_multislip_velbc.cfg
```

## Suggested exercises

1. Change the mesh to a quad mesh.

  * Create a journal file based on mesh_tri3.jou and use quads instead of tris.
  * Create new .cfg files (stepxx_yy.cfg) to use the quad mesh.

2. Change the material contrast across the fault.

3. Change the basis order and quadrature order for the solution field.

4. Change Step03 to include multiple ruptures at different time
   intervals, while still accounting for the accumulated slip deficit.
