# PyLith

[![MIT License](https://img.shields.io/badge/license-MIT-blue.svg)](https://github.com/geodynamics/pylith/blob/master/COPYING)
[![DOI](https://zenodo.org/badge/12650142.svg)](https://zenodo.org/badge/latestdoi/12650142)


## Build status
* Master: [![Build Status](https://travis-ci.org/geodynamics/pylith.svg?branch=master)](https://travis-ci.org/geodynamics/pylith)
* Maint: [![Build Status](https://travis-ci.org/geodynamics/pylith.svg?branch=maint)](https://travis-ci.org/geodynamics/pylith)
* Next: [![Build Status](https://travis-ci.org/geodynamics/pylith.svg?branch=next)](https://travis-ci.org/geodynamics/pylith)


PyLith is an open-source finite-element code for dynamic and
quasistatic simulations of crustal deformation, primarily earthquakes
and volcanoes.

* Main page: [https://geodynamics.org/cig/software/pylith](https://geodynamics.org/cig/software/pylith)
  * User Manual
  * Binary packages
  * Utility to build PyLith and all of its dependencies from source
* PyLith Wiki: [https://wiki.geodynamics.org/software:pylith:start](https://wiki.geodynamics.org/software:pylith:start)
  * Archive of online tutorials
  * Hints, tips, tricks, etc
  * PyLith development plan 
* Submit bug reports via https://github.com/geodynamics/pylith/issues
* Send all questions to: cig-short@geodynamics.org


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

See [CHANGES](CHANGES) for a complete list of changes for each release.

## Authors

* Brad Aagaard, Earthquake Science Center, USGS, USA
* Matthew Knepley, Computational and Applied Mathematics, Rice University, USA
* Charles Williams, Tectonophysics, GNS Science, New Zealand
