# Examples: 2-D reverse faults/gravity

This suite of examples demonstrates some basic concepts of using
PyLith with body forces as well as providing simple examples of
meshing and using faults with specified slip in 2-D.
Concepts common to all of the steps include:

* Dirichlet boundary conditions
* Isotropic, linear elasticity with a single material
* Output of the solution over the domain, boundaries, and materials
* Output of auxiliary information for boundary conditions and
  materials

## Meshing: Meshing a 2-D geometry including a main thrust and a splay fault

To do the meshing, Cubit/Trelis is required, and the journal files required are
`geometry.jou`, `createbc.jou`, `gradient.jou`, `mesh_tri3.jou` and
`mesh_quad4.jou`.

## Step01: Gravitational body forces with Dirichlet boundary conditions

Zero-displacement Dirichlet boundary conditions on the +x, -x, and
-y boundaries, and gravitational body forces applied. Features used in this
simulation include:

* Static simulation
* SimpleDB spatial database for specifying values for properties and
  boundary conditions
* Use of gravity_field to simulate gravitational body forces.

The simulation parameters are in the `pylithapp.cfg`,
`step01_gravity.cfg`, `step01_gravity_quad.cfg`  and 'step01_gravity_tri.cfg`
files. Note that we provide options for both triangular and quadrilateral cells.

To run the example (quad cells):
```
pylith step01_gravity.cfg step01_gravity_quad.cfg
```

To run the example (tri cells):
```
pylith step01_gravity.cfg step01_gravity_tri.cfg
```

## Step02: Gravitational body forces with reference stresses

Zero-displacement Dirichlet boundary conditions on the +x, -x, and
-y boundaries, gravitational body forces applied, and reference normal stresses
applied to balance the body forces. Features used in this simulation include:

* Static simulation
* SimpleDB spatial database for specifying values for properties and
  boundary conditions
* Use of gravity_field to simulate gravitational body forces.
* Use of reference_stress subfield to provide reference stresses.

The simulation parameters are in the `pylithapp.cfg`,
`step02_gravity_refstate.cfg`, `step02_gravity_refstate_quad.cfg`  and
'step02_gravity_refstate_tri.cfg` files. Note that we provide options for both
triangular and quadrilateral cells.

To run the example (quad cells):
```
pylith step02_gravity_refstate.cfg step02_gravity_refstate_quad.cfg
```

To run the example (tri cells):
```
pylith step02_gravity_refstate.cfg step02_gravity_refstate_tri.cfg
```

## Step03: Gravitational body forces and incompressible elasticity

The same as problem Step01, but in this case we use an incompressible
elasticity formulation combined with a large bulk modulus to simulate
incompressible conditions. Features used in this simulation include:

* Static simulation
* SimpleDB spatial database for specifying values for properties and
  boundary conditions
* Use of gravity_field to simulate gravitational body forces.
* Use of pylith.problems.SolnDispPres for a solution involving both
  displacements and pressure.
* Use of Dirichlet BC for the pressure field.

The simulation parameters are in the `pylithapp.cfg`,
`step03_gravity_incompressible.cfg`, `step03_gravity_incompressible_quad.cfg`
and 'step03_gravity_incompressible_tri.cfg` files. Note that we provide options
for both triangular and quadrilateral cells.

To run the example (quad cells):
```
pylith step03_gravity_incompressible.cfg step03_gravity_incompressible_quad.cfg
```

To run the example (tri cells):
```
pylith step03_gravity_incompressible.cfg step03_gravity_incompressible_tri.cfg
```

## Step04: Ground surface loading using Neumann BC

This problem imposes normal tractions on the ground surface, which could
simulate loading by water, ice, etc. Features used in this simulation include:

* Static simulation
* Neumann boundary conditions
* SimpleDB spatial database for specifying values for properties and
  boundary conditions.

The simulation parameters are in the `pylithapp.cfg`, `step04_neumann.cfg`,
`step04_neumann_quad.cfg` and 'step04_neumann_tri.cfg` files. Note that we
provide options for both triangular and quadrilateral cells.

To run the example (quad cells):
```
pylith step04_neumann.cfg step04_neumann_quad.cfg
```

To run the example (tri cells):
```
pylith step04_neumann.cfg step04_neumann_tri.cfg
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
