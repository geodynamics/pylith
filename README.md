# PyLith

[![DOI](https://www.zenodo.org/badge/DOI/10.5281/zenodo.7072811.svg)](https://doi.org/10.5281/zenodo.7072811)
[![MIT License](https://img.shields.io/badge/license-MIT-blue.svg)](https://github.com/geodynamics/pylith/blob/main/LICENSE.md)
[![Build Status](https://dev.azure.com/baagaard-usgs/pylith/_apis/build/status/geodynamics.pylith?branchName=main)](https://dev.azure.com/baagaard-usgs/pylith/_build/latest?definitionId=2&branchName=main)
[![codecov](https://codecov.io/gh/geodynamics/pylith/branch/master/graph/badge.svg?token=JiwLVB64EF)](https://codecov.io/gh/geodynamics/pylith)

## Description

PyLith is an open-source finite-element code for dynamic and
quasi-static simulations of crustal deformation, primarily earthquakes
and volcanoes.

* Main page: [https://geodynamics.org/resources/pylith](https://geodynamics.org/resources/pylith)
  * [User Guide](https://pylith.readthedocs.io/en/latest/user/index.html)
  * Binary packages
  * Utility to build PyLith and all of its dependencies from source
*
  [PyLith parameter viewer](https://geodynamics.github.io/pylith_parameters/) for viewing `.json` parameter files.
* PyLith Wiki: [https://wiki.geodynamics.org/software:pylith:start](https://wiki.geodynamics.org/software:pylith:start)
  * Archive of online tutorials
  * Hints, tips, tricks, etc
  * [PyLith development plan](https://github.com/geodynamics/pylith/wiki/Development-Plans)
* Submit bug reports via [GitHub issues](https://github.com/geodynamics/pylith/issues).
* Post all questions to the [PyLith category](https://community.geodynamics.org/c/pylith/29) on the [CIG Community Forum](https://community.geodynamics.org).


## Features

* Quasi-static (implicit) and dynamic (explicit) time-stepping
* Cell types include triangles, quadrilaterals, hexahedra, and tetrahedra
* Linear elastic, linear and generalized Maxwell viscoelastic, power-law viscoelastic, and Drucker-Prager elastoplastic materials
* Infinitesimal and small strain elasticity formulations
* Fault interfaces using cohesive cells
  * Prescribed slip with multiple, potentially overlapping earthquake ruptures and aseismic creep
  * Spontaneous slip with slip-weakening friction and Dieterich rate- and state-friction fault constitutive models
* Time-dependent Dirichlet (displacement/velocity) boundary conditions
* Time-dependent Neumann (traction) boundary conditions
* Time-dependent point forces
* Absorbing boundary conditions
* Gravitational body forces
* VTK and HDF5/Xdmf output of solution, fault information, and state variables
* Templates for adding your own bulk rheologies, fault constitutive models, and interfacing with a custom seismic velocity model.
* User-friendly computation of static 3-D Green's functions

## Release Notes

See [CHANGES](CHANGES.md) for a complete list of changes for each release.

## Authors

PyLith is continually being improved by a growing, collaborative, and inclusive community. It is primarily developed and maintained by:

* Brad Aagaard, Geologic Hazards Science Center, U.S. Geological Survey, USA
* Matthew Knepley, Computer Science and Engineering, University at Buffalo, USA
* Charles Williams, Tectonophysics, GNS Science, New Zealand

A complete and growing list of the authors that have contributed over the years can be found at [GitHub contributors](https://github.com/geodynamics/pylith/graphs/contributors).

Please see the User Guide for complete and [acknowledgement and citation](https://pylith.readthedocs.io/en/latest/intro/preface.html#citation) information.

