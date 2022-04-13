# Examples: 3-D box

This set of examples demonstrates some basic concepts of using
PyLith to solve the static or quasistatic elasticity equation in
3-D. Concepts common to all of the steps include:

* Mesh from CUBIT/Trelis
* Dirichlet boundary conditions
* Isotropic, linear elasticity with a single material
* Output of the solution over the domain, boundaries, and materials
* Output of auxiliary information for boundary conditions and
  materials

Faults are discussed in the 3d/strike-slip (planned) and 3d/subduction
example sets.

## Step01: Axial extension with Dirichlet boundary conditions

Axial extension with Dirichlet boundary conditions on the +x, -x, +y,
-y, and -z boundaries. Features used in this simulation include:

* Static simulation
* UniformDB spatial database for specifying values for properties and
  boundary conditions

The simulation parameters are in the `pylithapp.cfg` and
`step01_axialdisp.cfg` files.

To run the example:
```
pylith step01_axialdisp.cfg
```

## Step02:  Simple shear with Dirichlet boundary conditions

Simple shear with Dirichlet boundary conditions on the four lateral
boundaries and the bottom boundary. Features used in this simulation
include:

* Static simulation
* UniformDB and SimpleGridDB spatial database for specifying values
  for properties and boundary conditions

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

1. Change the mesh to the hex mesh.

  * Change the setting in the pylithapp.cfg file.
  * Override the current setting using the command line.
  
2.  Reduce the shear modulus while keeping the bulk modulus the same.

3. Change the basis order and quadrature order for the solution field.

4. Change Step01 to axial compression in the x and y directions.

5. Change Step03 to be a combination of shear and axial compression/extension.

6. In Step04 change the initial conditions so that only one of the components is equal to the solution.

7. In Step05 adjust the rate and/or time when the rate dependence starts.
