# Axial and Shear Deformation in a 3D Box

This suite of examples demonstrates some basic concepts of using PyLith to solve the static and quasistatic boundary elasticity equation in a 3D box with uniform material properties.
This example builds on the 2D version and incrementally adds complexity through a series of steps:

:Step 1: Axial extension with Dirichlet (displacement) boundary conditions.
:Step 2: Shear deformation with Dirichlet (displacement) boundary conditions.
:Step 3: Shear deformation with Dirichlet (displacement) and Neumann (traction) boundary conditions.
:Step 4: Same as Step 2 but with initial conditions equal to the analytical solution.
:Step 5: Shear deformation with time-dependent Dirichlet (displacement) and Neumann (traction) boundary conditions.

## Step 1: Axial extension with Dirichlet boundary conditions

Axial extension with Dirichlet boundary conditions on the +x, -x, +y,
-y, and -z boundaries.

To run the example:

```bash
pylith step01_axialdisp.cfg
```

## Step 2: Simple shear with Dirichlet boundary conditions

Simple shear with Dirichlet boundary conditions on the four lateral
boundaries and the bottom boundary.

To run the example:

```bash
pylith step02_sheardisp.cfg
```

## Step 3: Simple shear with Dirichlet and Neumann boundary conditions

The same problem as Step 2 but with Neumann (traction) boundary
conditions replacing the Dirichlet boundary conditions on the +y and
-y boundaries.

To run the example:

```bash
pylith step03_sheardisptract.cfg
```

## Step 4: Simple shear with Dirichlet boundary conditions and initial conditions

The same problem as Step02 but with the initial conditions set to the
solution.

To run the example:

```bash
pylith step04_sheardispic.cfg
```

## Step 5: Time-dependent shear with Dirichlet and Neumann boundary conditions

Similar to Step 3 but with time-dependent boundary conditions
consisting of an initial value and a rate of change that is added
starting at 1.0 year.  

To run the example:

```bash
pylith step05_sheardisptractrate.cfg
```
