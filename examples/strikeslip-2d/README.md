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

## Step 4: Variable Coseismic Slip

We use this example to illustrate prescribing slip that varies along the strike of the fault. This example also serves as a means to generate coseismic displacements at fake GPS stations. In Step 6 we will use the displacements at these stations along with static Green’s functions computed in Step 5 to invert for the slip on the fault.

To run the example:
```
pylith step04_varslip.cfg
```

## Step 5: Green's Functions

In this example we compute static Green’s functions for fault slip and use them in Step 6 to invert for fault slip.

To run the example:
```
pylith step05_greensfns.cfg
```

## Step 6: Least Squares Fault Slip Inversion

In this example we do a simple static slip inversion using least squares. We treat the displacements at the fake GPS stations in Step 4 as the “observations” and use the Green’s functions from Step 5 to invert for the fault slip that we prescribed in Step 4.

To run the example:
```
./invert_slip.py
```

## Step 7: Bayesian Fault Slip Inversion

Warning: This examples requires the CATMIP Bayesian inversion framework which is not yet publicly available. In this example we perform the same inversion as in Step 6, but replace the least squares inversion with the CATMIP Bayesian inversion framework. We demonstrate inverting for fault slip using both the original CATMIP algorithm (Step 7a) and the crossfade CATMIP algorithm (Step 7b).

To run the example using the original CATMIP algorithm:
```
mpiexec -n 2 catmip_pylith_staticslip step07a_catmip.in
```

To run the example using the crossfade CATMIP algorithm:
```
mpiexec -n 2 cfcatmip_pylith_staticslip step07b_cfcatmip.in
```
