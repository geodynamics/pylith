# Examples: 2-D box

This suite of examples demonstrates some basic concepts of using PyLith to solve the static or quasistatic elasticity equation in 2-D. Faults are discussed in the 2d/strike-slip and 2d/reverse example suites.



## Step01: Axial extension with Dirichlet boundary conditions

Axial extension with Dirichlet boundary conditions on the +x, -x, and -y boundaries.
The simulation parameters are in the `pylithapp.cfg` and `step01_axialdisp.cfg` files.

To run the example:
```
pylith step01_axialdisp.cfg
```

## Step02:  Simple shear with Dirichlet boundary conditions

Simple shear with Dirichlet boundary conditions on all four boundaries.
The simulation parameters are in the `pylithapp.cfg` and `step01_sheardisp.cfg` files.

To run the example:
```
pylith step02_sheardisp.cfg
```

## Step03: Simple shear with Dirichlet and Neumann boundary conditions

The same problem as Step02 but with Neumann (traction) boundary conditions replacing the Dirichlet boundary conditions on the +x and +y boundaries.  The simulation parameters are in the `pylithapp.cfg` and `step01_sheardisptract.cfg` files.

To run the example:
```
pylith step03_sheardisptract.cfg
```

## Step04: Simple shear with Dirichlet and Neumann boundary conditions and initial conditions

The same problem as Step03 but with the initial conditions set to the solution. The simulation parameters are in the `pylithapp.cfg` and `step01_sheardisptractic.cfg` files.

To run the example:
```
pylith step04_sheardisptractic.cfg
```

## Step05: Time-dependent shear with Dirichlet and Neumann boundary conditions

Similar to Step03 but with time-dependent boundary conditions consisting of an initial value and a rate of change that is added starting at 1.0 year.  The simulation parameters are in the `pylithapp.cfg` and `step01_sheardisptractrate.cfg` files.

To run the example:
```
pylith step05_sheardisptract.cfg
```

## Suggested exercises

1. Change the mesh to the tri mesh.

  * Change the setting in the pylithapp.cfg file.
  * Override the current setting using the command line.
  
2.  Reduce the shear modulus while keeping the bulk modulus the same.

3. Change the basis order and quadrature order for the solution field.

4. Change Step01 to axial compression in the +y direction.

5. Change Step03 to be a combination of shear and axial compression/extension.

6. In Step04 change the initial conditions so that only one of the components is equal to the solution.

7. In Step05 adjust the rate and/or time when the rate dependence starts.
