# Development Plan

Future implementation of features is guided by several target applications, including

* Earthquake cycle modeling with quasistatic simulation of interseismic deformation and dynamic simulation of coseismic deformation.
* Inversion of geodetic data for slow slip events, fault creep, and long-term fault slip rates.
* Quasistatic and dynamic modeling of fluids and faulting.

## Version 3.0.1 (June 2022)

Finish updating examples and documentation.

## Version 3.1.0 (October 2022)

* Parallel mesh loading ![expert](images/expert.png) [15%]
* Diagonal Jacobian for explicit part of IMEX formulation ![expert](images/expert.png)[0%]
* Improve robustness of HDF5 output by opening/closing at each time step ![easy](images/easy.png)[0%]
* Additional minor cleanup of code internals to improve maintainability.

## Version 3.2 (December 2022)

* Spontaneous rupture for quasistatic and dynamic simulations ![expert](images/expert.png) [20%]
* Reimplementation of small strain formulation for elasticity ![difficult](images/difficult.png) [20%]
* Reimplementation of Drucker-Prager elastoplastic bulk rheology ![intermediate](images/intermediate.png) [0%]
* Add support for GeoModelGrids implementation of spatial databases for 3D seismic velocity models.
* Line/point fluid sources in poroelasticity ![expert](images/expert.png) [20%]
* Integration with libCEED for fast high order residual evaluation ![expert](images/expert.png)\
  Contribution led by Jed Brown.
* Add ability to output residual field during nonlinear solve for debugging ![easy](images/easy.png) [0%]

## Version 3.3 (June 2023)

* Coupling of problems with compatible meshes ![difficult](images/difficult.png) [10%]\
    Implement "injectors" for solution and state variables.
* Migrate examples to Jupyter notebooks ![intermediate](images/intermediate.png)
* Update to current version of Pyre ![difficult](images/difficult.png)
* Moment tensor point sources  [5%]\
  Moment tensor point sources provide a mesh independent deformation source that is better suited for Green's function calculations than slip on a fault surface via cohesive cells.
* Adaptive mesh refinement ![expert](images/expert.png)

## Features for Future Releases

* Consolidate HDF5 output into a single file ![difficult](images/difficult.png)
* Elasticity with self-gravitation ![intermediate](images/intermediate.png)
* Drucker-Prager bulk rheology with relaxation to yield surface ![intermediate](images/intermediate.png) 
* Drucker-Prager bulk rheology with strain hardening/softening  ![intermediate](images/intermediate.png)
* Adjoint for data assimilation ![difficult](images/difficult.png)
* Fault with both prescribed slip and spontaneous rupture ![difficult](images/difficult.png)\
  Use fault constitutive model to control slip on fault except during episodes of prescribed slip. Need some way to describe when to turn on/off prescribed slip.

