# Development Plan

Future implementation of features is guided by several target applications, including

* Earthquake cycle modeling with quasi-static simulation of interseismic deformation and dynamic simulation of coseismic deformation.
* Inversion of geodetic data for slow slip events, fault creep, and long-term fault slip rates.
* Quasistatic and dynamic modeling of fluids and faulting.

## Version 3.0.4 (September 2023)

Updates for examples and documentation along with bugfixes.

* Analytical function spatial database ![easy](images/easy.png)[100%]
* Use `pod` initial guess to improve convergence ![easy](images/easy.png) [100%]
* Dynamic prescribed slip with diagonal Jacobian for explicit part of IMEX formulation ![expert](images/expert.png)[75%]
* Output of fault tractions ![expert](images/expert.png) [75%]
* Update coordinates with solution ![intermediate](images/intermediate.png) [50%]
* Convert from CppUnit to Catch2 ![easy](images/easy.png) [40%]
* Better preconditioners ![expert](images/expert.png) [25%]
  * elasticity with fault
  * incompressible elasticity
  * poroelasticity
* Finish updating `examples/subduction-3d` ![intermediate](images/intermediate.png) [20%]
* Add `examples/barwaves-2d` ![expert](images/expert.png) [15%]
* Parallel mesh loading ![expert](images/expert.png) [15%]
* Elasticity with self-gravitation ![intermediate](images/intermediate.png) [0%]
* Add 2D and 3D examples for crustal faults with complex fault geometry ![easy](images/easy.png) [0%]

## Version 3.1.0 (December 2023)

* Spontaneous rupture for quasistatic simulations ![expert](images/expert.png) [20%]
* Add support for GeoModelGrids implementation of spatial databases for 3D seismic velocity models. ![intermediate](images/intermediate.png) [0%]
* Improve robustness of HDF5 output by opening/closing at each time step ![easy](images/easy.png)[0%]
* Dirichlet boundary conditions with constraints on normal and tangential components. ![difficult](images/difficult.png) [0%]
* Additional minor cleanup of code internals to improve maintainability.
  * Refactor auxiliary field and derived field output (initial/timestep/final)
  * Output of fault rupture auxiliary fields
  * VTK output (vtk -> vtu)
* Reimplementation of Drucker-Prager elastoplastic bulk rheology ![intermediate](images/intermediate.png) [0%]
* Convert from SWIG to pybind11 ![intermediate](images/intermediate.png) [0%]

## Version 4.0 (June 2024)

* Spontaneous rupture for dynamic simulations ![expert](images/expert.png) [10%]
* Migrate examples to Jupyter notebooks ![intermediate](images/intermediate.png)
* Update to current version of Pyre ![difficult](images/difficult.png)
* More flexible specification of time-dependent boundary conditions. ![difficult](images/difficult.png) [0%]
* Integration with libCEED for fast high order residual evaluation ![expert](images/expert.png)\
  Contribution led by Jed Brown.
* Add ability to output residual field during nonlinear solve for debugging ![easy](images/easy.png) [0%]

## Features for Future Releases

* Coupling of problems with compatible meshes ![difficult](images/difficult.png) [10%]\
    Implement "injectors" for solution and state variables.
* Reimplementation of small strain formulation for elasticity ![difficult](images/difficult.png) [20%]
* Moment tensor point sources  [5%]\
  Moment tensor point sources provide a mesh independent deformation source that is better suited for Green's function calculations than slip on a fault surface via cohesive cells.
* Adaptive mesh refinement ![expert](images/expert.png)
* Line/point fluid sources in poroelasticity ![expert](images/expert.png) [20%]
* Consolidate HDF5 output into a single file ![difficult](images/difficult.png)
* Drucker-Prager bulk rheology with relaxation to yield surface ![intermediate](images/intermediate.png) 
* Drucker-Prager bulk rheology with strain hardening/softening  ![intermediate](images/intermediate.png)
* Adjoint for data assimilation ![difficult](images/difficult.png)
* Fault with both prescribed slip and spontaneous rupture ![difficult](images/difficult.png)\
  Use fault constitutive model to control slip on fault except during episodes of prescribed slip. Need some way to describe when to turn on/off prescribed slip.
