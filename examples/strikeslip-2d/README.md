# Examples: 2D Strike-Slip Fault

This suite of examples demonstrates some basic concepts of using
PyLith with a single through-going strike-slip fault and a simple
rectangular domain.

## Meshing

We provide mesh files generated using Gmsh and Cubit.
We also include a Python script for generating the finite-element mesh with
triangular cells using Gmsh and Journal files for generating the
finite-element mesh with triangular or quadrilateral cell using Cubit.

## Step 1: Static Coseismic Slip

This example involves a static simulation that solves for the deformation from prescribed coseismic slip on the fault. We specify 2 meters of right-lateral slip.

To run the example:
```
pylith step01_slip.cfg
```

## Step 2: Single Earthquake Rupture and Velocity Boundary Conditions

This example involves a quasistatic simulation that solves for the deformation from velocity boundary conditions and prescribed coseismic slip on the fault.
We let strain accumulate due to the motion of the boundaries and then release the strain by prescribing 2 meters of right-lateral slip at t=100 years.

To run the example:
```
pylith step02_slip_velbc.cfg
```

## Step 3: Multiple Earthquake Ruptures and Velocity Boundary Conditions

This example involves a quasistatic simulation that solves for the deformation from velocity boundary conditions and multiple earthquake ruptures on the fault.
The velocity boundary conditions match those in Step 2.
We prescribe the first earthquake rupture to occur at 100 years with 1 meter of right-lateral slip and the second earthquake rupture to occur at 200 years with 3 meters of right-lateral slip.

To run the example:
```
pylith step03_multislip_velbc.cfg
```
