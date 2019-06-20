# Examples: 2-D Reverse Fault and Gravitational Body Forces

This suite of examples demonstrates some basic concepts of using
PyLith with body forces as well as providing simple examples of
meshing and using faults with prescribed slip in 2-D.
Concepts common to all of the steps include:

* Dirichlet boundary conditions
* Isotropic, linear elasticity with multiple materials but identical
  properties
* Output of the solution over the domain, boundaries, and materials
* Output of auxiliary information for boundary conditions and
  materials

## Meshing: Meshing a 2-D geometry including a main thrust and a splay fault

We create the mesh using CUBIT/Trelis; however the resulting mesh is
provided, so you can skip creating the mesh if you do not have
CUBIT/Trelis. The associated journal files are: `geometry.jou`,
`createbc.jou`, `gradient.jou`, `mesh_tri.jou` (triangule cells) and
`mesh_quad.jou` (quadrilateral mesh).

## Step01: Gravitational body forces with Dirichlet boundary conditions

Zero-displacement Dirichlet boundary conditions on the +x, -x, and -y
boundaries, and gravitational body forces applied. Features used in
this simulation include:

* Static simulation
* SimpleDB and ZeroDB spatial database for specifying values for
  properties and boundary conditions
* Use of gravity_field to simulate gravitational body forces

The simulation parameters are in the `pylithapp.cfg` and
`step01_gravity.cfg`

To run the example:
```
pylith step01_gravity.cfg
```

## Step02: Gravitational body forces with reference stresses

Zero-displacement Dirichlet boundary conditions on the +x, -x, and
-y boundaries, gravitational body forces applied, and reference normal stresses
applied to balance the body forces. Features used in this simulation include:

* Static simulation
* SimpleDB and ZeroDB spatial database for specifying values for
  properties and boundary conditions
* Use of gravity_field to simulate gravitational body forces.
* Use of reference_stress subfield to provide reference stresses.

The simulation parameters are in the `pylithapp.cfg` and `step02_gravity_refstate.cfg` files.

To run the example:
```
pylith step02_gravity_refstate.cfg
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
* Additional solver parameters to solve an incompressible elasticity problem.

The simulation parameters are in the `pylithapp.cfg` and `step03_gravity_incompressible.cfg` files.

To run the example:
```
pylith step03_gravity_incompressible.cfg
```

## Step04: Ground surface loading using Neumann BC

This problem imposes normal tractions on the ground surface, which could
simulate loading by water, ice, etc. Features used in this simulation include:

* Static simulation
* Neumann boundary conditions
* SimpleDB and ZeroDB spatial database for specifying values for
  properties and boundary conditions.

The simulation parameters are in the `pylithapp.cfg` and `step04_surfload.cfg` files.

To run the example:
```
pylith step04_surfload.cfg
```

## Step05: Fault slip on a single fault

This problem has the same zero-displacement BC used in Step 01-03, but
includes earthquake rupture of the main thrust fault. Features used in
this simulation include:

* Static simulation
* UniformDB, ZeroDB, and SimpleDB spatial database for specifying
  values for properties, faults, and boundary conditions
* Use of pylith.problems.SolnDispLagrange for problems involving
  faults
* Additional solver parameters to solve a problem involving faults.

The simulation parameters are in the `pylithapp.cfg`, `step05_onefault.cfg`
and `solver_faults.cfg` files. Note that the additional solver settings in
'solver_faults.cfg' are required for all of the faulting examples
(step05-step08).

To run the example:
```
pylith step05_onefault.cfg solver_faults.cfg
```

## Step06: Fault slip on main fault and splay fault

This problem is similar to Step05, but includes earthquake rupture on both
the main thrust fault and a splay fault. Features used in this simulation
include:

* Static simulation
* UniformDB, ZeroDB, and SimpleDB spatial database for specifying
  values for properties, faults, and boundary conditions
* Use of pylith.problems.SolnDispLagrange for problems involving
  faults
* Slip on multiple faults, with the main thrust considered the throughgoing
  fault
* Additional solver parameters to solve a problem involving faults.

The simulation parameters are in the `pylithapp.cfg`, `step06_twofault.cfg`
and `solver_faults.cfg` files. Note that the additional solver settings in
'solver_faults.cfg' are required for all of the faulting examples
(step05-step08).

To run the example:
```
pylith step06_twofault.cfg solver_faults.cfg
```

## Step07: Fault slip on main fault and splay fault with Maxwell viscoelasticity

**Note**: This example does not work in v3.0.0beta3. It will be
working in the v3.0.0 release.

This problem is similar to Step06, but the slab material is considered to be
viscoelastic, and the simulation runs for 100 years. Features used in this
simulation include:

* Quasi-static simulation
* UniformDB, ZeroDB, and SimpleDB spatial database for specifying
  values for properties, faults, and boundary conditions
* Use of pylith.problems.SolnDispLagrange for problems involving
  faults
* Additional solver parameters to solve a problem involving faults.
* Slip on multiple faults, with the main thrust considered the throughgoing
  fault
* Use of a Maxwell viscoelastic material
* Time-dependent problem running for 100 years

The simulation parameters are in the `pylithapp.cfg`,
`step07_twofault_maxwell.cfg` and `solver_faults.cfg` files. Note that the
additional solver settings in 'solver_faults.cfg' are required for all of
the faulting examples (step05-step08).

To run the example:
```
pylith step07_twofault_maxwell.cfg solver_faults.cfg
```

## Step08: Fault slip on main fault and splay fault with power-law viscoelasticity

**Note**: This example does not work in v3.0.0beta3. It will be
working in the v3.0.0 release.

This problem is similar to Step07, but the slab material is considered to be
power-law viscoelastic, and the simulation runs for 100 years. Features used
in this simulation include:

* Quasi-static simulation
* UniformDB, ZeroDB, and SimpleDB spatial database for specifying
  values for properties, faults, and boundary conditions
* Use of pylith.problems.SolnDispLagrange for problems involving
  faults
* Slip on multiple faults, with the main thrust considered the throughgoing
  fault
* Additional solver parameters to solve a problem involving faults.
* Use of a power-law (nonlinear) viscoelastic material
* Time-dependent problem running for 100 years
* Required use of the nonlinear solution type

The simulation parameters are in the `pylithapp.cfg`,
`step08_twofault_powerlaw.cfg` and `solver_faults.cfg` files. Note that the
additional solver settings in 'solver_faults.cfg' are required for all of
the faulting examples (step05-step08).

To run the example:
```
pylith step08_twofault_powerlaw.cfg solver_faults.cfg
```

## Suggested exercises

1. Change the mesh to the quad mesh.

  * Change the setting in the pylithapp.cfg file.
  * Override the current setting using the command line.
  
2. Use different elastic properties/density for the elastic problems (Step01-06).

3. Change the basis order and quadrature order for the solution field.

4. Change the shape of the applied surface load in Step04.

5. Try using variable fault slip for the faulting problems (Step05-08).

6. Introduce multiple ruptures on one of the faults for the faulting problems
   (Step05-08).
