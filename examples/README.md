All:
  HDF5 output

2d/box
  tri mesh
  quad mesh (default)
  Material: single, isotropic linear elasticity; uniform properties
  Output: domain, material, bc
  basis order (default: 1); quad order (default: 1)

  Step01: axial extension w/Dirichlet BC
    Dirichlet BC (UniformDB)
    IsotropicLinearElasticity (UniformDB)
  Step02: simple shear w/Dirichlet BC
    Dirichlet BC (SimpleDB)
    IsotropicLinearElasticity (UniformDB)
  Step03: simple shear w/Dirichlet + Neumann BC
    Dirichlet BC (SimpleDB)
    Neumann BC (UniformDB)
    IsotropicLinearElasticity (SimpleDB)
  Step04: simple shear w/Dirichlet + Neumann BC + initial conditions
    Dirichlet BC (SimpleDB)
    Neumann BC (UniformDB)
    IsotropicLinearElasticity (SimpleDB)
    IC: SimpleDB
  Step03: time-dependent shear w/Dirichlet and Neumann BC
    Dirichlet BC (UniformDB)
    Neumann BC (UniformDB)
    IsotropicLinearElasticity (SimpleDB)
  Exercises
    Change to tri mesh
    Change material properties
1    Change BC to give axial compression in the +y direction
    Change basis and quadrature order

3d/box
  tet mesh (default)
  hex mesh
  two materials (upper & lower crust)

  Step01: axial extension w/Dirichlet BC
    Dirichlet BC (UniformDB)
    IsotropicLinearElasticity (2x UniformDB)
  Step02: simple shear w/Dirichlet BC
    Dirichlet BC (SimpleDB)
    IsotropicLinearElasticity (1x SimpleDB)
  Step03: simple shear w/Dirichlet + Neumann BC
    Dirichlet BC (SimpleDB)
    Neumann BC (UniformDB)
    IsotropicLinearElasticity (1x SimpleDB)
  Step03: time-dependent shear w/Dirichlet and Neumenn BC
    Dirichlet BC (UniformDB)
    Neumann BC (UniformDB)
    IsotropicLinearElasticity (1x SimpleDB)
  Exercises
    Change to hex mesh
    Change material properties
    Change BC to give axial compression in the +y direction
    Change basis and quadrature order


2d/strike-slip w/throughgoing fault
  tri mesh
  quad mesh (default)
  two materials (25% contrast)

  Step01: static w/prescribed slip; fixed boundaries
  Step02: quasistatic: boundaries move; fault catches up w/prescribed slip
  Step03: quasistatic: multiple earthquake cycles w/presscribed slip

  Exercises
    Change to tri mesh

2d/reverse w/buried fault + splay fault
  tri mesh (default)
  quad mesh
  three materials: upper crust, wedge, lower crust

  Step01: gravitational body forces (no reference state)
  Step02: gravitational body forces (reference state)
    Reference state (SimpleDB)
  Step03: gravitational body forces + incompressible elasticity
  Step04: gravitational body forces + incompressible elasticity
    Add initial condition to Step03 (SimpleDB)
  Step05: distributed surface load
    Neumann BC (SimpleGridDB)
  Step06: coseismic slip (main fault)
  Step07: coseismic slip (main + splay fault)
  Step08: coseismic slip + viscoelastic relaxation w/linear Maxwell
  Step09: coseismic slip + viscoelastic relaxation w/powerlaw

  Exercises
    Change to quad mesh
    Compare basis order 1 and 2 for Steps 1 and 2
    Change material properties

3d/strike-slip w/throughgoing fault :LATER:?
  match 2d/strike-slip w/throughgoing fault

2d/subduction :UPDATE:

2d/magmachamber :LATER:

3d/subduction :UPDATE:

3d/strike-slip :LATER:

debugging :UPDATE:
  quadrature orders don't match
  typos in component names

meshing [no changes]

Concepts:
  * Output
    + domain
      - solution
    + boundary
      - solution
    + fault
    + material
      - auxiliary subfields
      - derived field
      - solution
    + points
      - solution
    + Dirichlet BC
      - auxiliary field
    + Neumann BC
      - auxiliary field
    + Add field to HDF5 and update Xdmf
    + Projection of field to basis order 0/1
  * Spatial databases
    + Generate spatial database using Python
    + SimpleGridDB
    + UniformDB
  * Simulation
    + Neumann BC
      - side loads
      - surface loads
    + Dirichlet BC
      - roller BC
      - constant rate BC
    + Initial conditions
    + Elasticity
      - Referense state
      - Body forces
      - Gravitational body forces
      - Rheologies
        * IsotropicLinearElasticity
	* IsotropicLinearMaxwell
	* IsotropicPowerlaw
	* IsotropicDruckerPrager
    + Incompressible elasticity
      - gravitational body forces
      - initial conditions
    + Fault
      + Prescribed slip step function, single rupture
      + Prescribed slip step function + creep, multiple ruptures
  * Order of basis functions
    - gravity w/1st order vs 2nd order
