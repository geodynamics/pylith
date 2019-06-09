# Troubleshooting: 2-D box

This suite of examples demonstrates how to use a variety of techniques
to troubleshoot errors in running PyLith simulations.  The input files
contain several common mistakes and a few uncommon mistakes.

## Step01: Axial extension with Dirichlet boundary conditions

Axial extension with Dirichlet boundary conditions on the +x, -x, and
-y boundaries. Features used in this simulation include:

The simulation parameters are in the `pylithapp.cfg` and
`step01_axialdisp.cfg` files.

To run the example:
```
pylith step01_axialdisp.cfg
```

## Step02:  Simple shear with Dirichlet boundary conditions

Simple shear with Dirichlet boundary conditions on all four
boundaries. Features used in this simulation include:

The simulation parameters are in the `pylithapp.cfg` and
`step02_sheardisp.cfg` files.

To run the example:
```
pylith step02_sheardisp.cfg
```

## Step03: Simple shear with Dirichlet and Neumann boundary conditions

The same problem as Step02 but with Neumann (traction) boundary
conditions replacing the Dirichlet boundary conditions on the +x and
+y boundaries.  Features used in this simulation include:

The simulation parameters are in the `pylithapp.cfg` and
`step03_sheardisptract.cfg` files.

To run the example:
```
pylith step03_sheardisptract.cfg
```

