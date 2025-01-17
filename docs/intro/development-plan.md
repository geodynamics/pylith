# Development Plan

Future implementation of features is guided by several target applications, including

* Earthquake cycle modeling with quasi-static simulation of interseismic deformation and dynamic simulation of coseismic deformation.
* Inversion of geodetic data for slow slip events, fault creep, and long-term fault slip rates.
* Quasistatic and dynamic modeling of fluids and faulting.

We also strive to find the correct balance between adding new features and improving the code internals, performance, and maintainability.
Short-term priorities focus on reimplementing features in version 2.2 that have not yet been migrated to the current formulations for governing equations and discretization.

:::{figure-md} fig:development:plan
<img src="figs/development-1.*" alt="Diagram of development priorities."  width="100%"/>

PyLith development priorities by category.
Information to the right of each feature indicates difficulty (colored icon), the relative amount of effort (lower left), and an estimate of how much is completed (lower right percentage).
The arrows indicate the order in which some features must be implemented.
:::

## Anticipated Releases

:::{note}
Because we strictly follow the [semantic versioning guidelines](https://semver.org/), a minor release may get promoted to a major releases if we make changes to the public API (parameters).
This can happen if realize that we should modify the parameters to improve maintainability or prepare for future changes.
:::

## Version 4.3 (April 2025)

* Specify boundary conditions using faces instead of vertices ![expert](images/intermediate.png) [75%]
* Parallel mesh loading ![expert](images/expert.png) [85%]

## Version 4.4 (June 2025)

* Dynamic prescribed slip with diagonal Jacobian for explicit part of IMEX formulation ![expert](images/expert.png) [80%]
* Add `examples/barwaves-2d` ![expert](images/expert.png) [25%]
* Better preconditioner for incompressible elasticity ![expert](images/expert.png) [80%]

## Version 4.5 (October 2025)

* Spontaneous rupture for quasistatic and dynamic simulations ![expert](images/expert.png) [30%]


## Version 5.0 (June 2026)

* Convert from SWIG to pybind11 ![intermediate](images/intermediate.png) [0%]
* Add support for GeoModelGrids implementation of spatial databases for 3D seismic velocity models. ![intermediate](images/intermediate.png) [0%]
* Update to current version of Pyre ![difficult](images/difficult.png)
* Migrate examples to Jupyter notebooks ![intermediate](images/intermediate.png)
* Update VTK output to use `vtu` files rather than legacy `vtk` files ![easy](images/easy.png) [0%]
* Improve robustness of HDF5 output by opening/closing at each time step ![easy](images/easy.png)[0%]

## Version 6.0 (TBD)

* Output of fault rupture auxiliary subfields ![intermediate](images/intermediate.png) [0%]
* Improve creation of auxiliary, diagnostic, and derived fields.
* Reimplement Drucker-Prager elastoplastic bulk rheology ![intermediate](images/intermediate.png) [0%]

## Version 7.0.0 (TBD)

* Dirichlet boundary conditions with constraints on normal and tangential components. ![difficult](images/difficult.png) [0%]
* Integration with libCEED for fast high order residual evaluation ![expert](images/expert.png)\
* Add ability to output residual field during nonlinear solve for debugging ![easy](images/easy.png) [0%]

## Features for Future Releases

These features are not yet assigned to a release.
They may be added to the planned releases based upon further discussion and contributions from the community.

* Adaptive mesh refinement ![expert](images/expert.png)
* Faults with poroelastic properties ![expert]{images/expert.png} [80%]
* More flexible specification of time-dependent boundary conditions. ![difficult](images/difficult.png) [0%]
* Elasticity with self-gravitation ![intermediate](images/intermediate.png) [0%]
* Coupling of problems with compatible meshes ![difficult](images/difficult.png) [10%]\
    Implement "injectors" for solution and state variables.
* Reimplementation of small strain formulation for elasticity ![difficult](images/difficult.png) [20%]
* Moment tensor point sources  [60%]\
  Moment tensor point sources provide a mesh independent deformation source that is better suited for Green's function calculations than slip on a fault surface via cohesive cells.
* Line/point fluid sources in poroelasticity ![expert](images/expert.png) [20%]
* Consolidate HDF5 output into a single file ![difficult](images/difficult.png)
* Drucker-Prager bulk rheology with relaxation to yield surface ![intermediate](images/intermediate.png)
* Drucker-Prager bulk rheology with strain hardening/softening  ![intermediate](images/intermediate.png)
* Adjoint for data assimilation ![difficult](images/difficult.png)
* Fault with both prescribed slip and spontaneous rupture ![difficult](images/difficult.png)\
  Use fault constitutive model to control slip on fault except during episodes of prescribed slip. Need some way to describe when to turn on/off prescribed slip.
