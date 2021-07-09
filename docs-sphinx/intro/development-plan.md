# Development Plan

## Version 3.0.0dev

List of new features currently available in the `main` branch.

* Multiphysics
  * Elasticity for linear isotropic materials and linear Maxwell and generalized Maxwell models
  * Incompressible elasticity via a pressure field
  * Prescribed slip for quasistatic simulations
* Higher order basis functions\
    Allow user to select order of basis functions independent of the mesh (which defines the geometry). This permits higher resolution for a given mesh.
* Switch to using PETSc time-stepping (TS) algorithms\
  Replace simple Python-based time-stepping implementations with PETSc time-stepping algorithms that provide support for higher order discretization in time and real adaptive time stepping.
* Modular approach for initial conditions
* Output of subfields with user-defined basis order
* Simulation metadata with command line utility for searching metadata
* Convert to Python 3

## Version 3.0.0

* Multiphysics
  * IMEX/DAE implementation for prescribed slip for dynamic simulations. ![expert](images/expert.png) [80%]
  * Finish reimplementation of powerlaw viscoelastic bulk rheology ![intermediate](images/intermediate.png) [90%]
  * Poroelasticity with bulk isotropic, linear elasticity bulk rheology ![expert](images/expert.png) [85%]\
  Finish documentation (discussion of formulation in manual along with examples)
* Static Green's functions ![intermediate](images/intermediate.png) [75%]
* Ability to run multiple, independent problems in a single simulation [75%]
* Testing via Method of Manufactured Solutions ![intermediate](images/intermediate.png) [90%]

## Version 3.1 (Fall 2021)

* Spontaneous rupture for quasistatic and dynamic simulations ![expert](images/expert.png) [20%]
* Coupling of problems with compatible meshes ![difficult](images/difficult.png) [10%]\
    Implement `Injectors` for solution and state variables.
* Reimplementation of small strain formulation for elasticity ![difficult](images/difficult.png) [20%]
* Reimplementation of Drucker-Prager elastoplastic bulk rheology ![intermediate](images/intermediate.png) [0%]
* Line/point fluid sources in poroelasticity ![expert](images/expert.png) [20%]
* Parallel mesh loading ![expert](images/expert.png) [5%]\
  Requires creating cohesive cells in parallel.
* PETSc mesh importer, including support for Gmsh ![easy](images/easy.png) [30%]
* Integration with libCEED for fast high order residual evaluation ![expert](images/expert.png)\
  Contribution led by Jed Brown.
* Add ability to output residual field during nonlinear solve for debugging ![easy](images/easy.png) [0%]
* Finish migrating documentation from LaTeX to Sphinx+MyST

## Version 3.2 (Spring 2022)

* Migrate examples to Jupyter notebooks
* Moment tensor point sources  [5%]\
  Moment tensor point sources provide a mesh independent deformation source that is better suited for Green's function calculations than slip on a fault surface via cohesive cells.
* Add support for GeoModelGrids implementation of spatial databases for 3D seismic velocity models.

## Features for Future Releases

* Consolidate HDF5 output into a single file ![difficult](images/difficult.png)
* Elasticity with self-gravitation ![intermediate](images/intermediate.png)
* Drucker-Prager bulk rheology with relaxation to yield surface ![intermediate](images/intermediate.png) 
* Drucker-Prager bulk rheology with strain hardening/softening  ![intermediate](images/intermediate.png)
* Adaptive mesh refinement ![expert](images/expert.png)
* Adjoint for data assimilation ![difficult](images/difficult.png)
* Fault with both prescribed slip and spontaneous rupture ![difficult](images/difficult.png)\
  Use fault constitutive model to control slip on fault except during episodes of prescribed slip. Need some way to describe when to turn on/off prescribed slip.

