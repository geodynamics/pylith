# Examples: 2D box

This suite of examples demonstrates some basic concepts of using PyLith to solve the static and quasistatic boundary elasticity equation in a 2D box with uniform material properties.
This example incrementally adds complexity through a series of steps:

:Step 1: Axial deformation with Dirichlet (displacement) boundary conditions.
:Step 2: Shear deformation with Dirichlet (displacement) boundary conditions.
:Step 3: Shear deformation with Dirichlet (displacement) and Neumann (traction) boundary conditions.
:Step 4: Same as Step 3 but with initial conditions equal to the analytical solution.
:Step 5: Shear deformation with time-dependent Dirichlet (displacement) and Neumann (traction) boundary conditions.

## Step 1: Axial extension with Dirichlet boundary conditions

Axial extension with Dirichlet boundary conditions on the +x, -x, and
-y boundaries.

To run the example:
```
pylith step01_axialdisp.cfg
```

## Step 2:  Simple shear with Dirichlet boundary conditions

Simple shear with Dirichlet boundary conditions on all four boundaries.

To run the example:
```
pylith step02_sheardisp.cfg
```

## Step03: Simple shear with Dirichlet and Neumann boundary conditions

The same problem as Step02 but with Neumann (traction) boundary
conditions replacing the Dirichlet boundary conditions on the +x and
+y boundaries.  Features used in this simulation include:

* Static simulation
* Neumann boundary conditions
* UniformDB and SimpleDB spatial database for specifying values
  for properties and boundary conditions

The simulation parameters are in the `pylithapp.cfg` and
`step03_sheardisptract.cfg` files.

To run the example:
```
pylith step03_sheardisptract.cfg
```

## Step04: Simple shear with Dirichlet boundary conditions and initial conditions

The same problem as Step02 but with the initial conditions set to the
solution. Features used in this simulation include:

* Static simulation
* Neumann boundary conditions
* Initial conditions over the domain
* UniformDB, SimpleGridDB, and SimpleDB spatial database for
  specifying values for properties, boundary conditions, and initial
  conditions.

The simulation parameters are in the `pylithapp.cfg` and
`step04_sheardispic.cfg` files.

To run the example:
```
pylith step04_sheardispic.cfg
```

## Step05: Time-dependent shear with Dirichlet and Neumann boundary conditions

Similar to Step03 but with time-dependent boundary conditions
consisting of an initial value and a rate of change that is added
starting at 1.0 year.  Features used in this simulation include:

* Quasistatic simulation
* Time-dependent Dirichlet boundary conditions
* Time-dependent Neumann boundary conditions
* UniformDB and SimpleDB spatial database for specifying values for
  properties and boundary conditions.

The simulation parameters are in the `pylithapp.cfg` and
`step05_sheardisptractrate.cfg` files.

To run the example:
```
pylith step05_sheardisptractrate.cfg
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
