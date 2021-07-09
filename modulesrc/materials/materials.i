// ======================================================================
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University at Buffalo
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2021 University of California, Davis
//
// See LICENSE.md for license information.
//
// ======================================================================
//

// SWIG interface
%module materials

// Header files for module C++ code
%{
#include "pylith/materials/Material.hh"
#include "pylith/materials/Elasticity.hh"
#include "pylith/materials/RheologyElasticity.hh"
#include "pylith/materials/IsotropicLinearElasticity.hh"
#include "pylith/materials/IsotropicLinearMaxwell.hh"
#include "pylith/materials/IsotropicLinearGenMaxwell.hh"
#include "pylith/materials/IsotropicPowerLaw.hh"
#include "pylith/materials/IncompressibleElasticity.hh"
#include "pylith/materials/RheologyIncompressibleElasticity.hh"
#include "pylith/materials/IsotropicLinearIncompElasticity.hh"
#include "pylith/materials/Poroelasticity.hh"
#include "pylith/materials/RheologyPoroelasticity.hh"
#include "pylith/materials/IsotropicLinearPoroelasticity.hh"

#include "pylith/utils/arrayfwd.hh"
%}

%include "exception.i"
%exception {
  try {
    $action
  } catch (const std::exception& err) {
    SWIG_exception(SWIG_RuntimeError, err.what());
  } // try/catch
 } // exception

%include "typemaps.i"
%include "../include/scalartypemaps.i"

// Numpy interface stuff
%{
#define SWIG_FILE_WITH_INIT
%}
%include "../include/numpy.i"
%init %{
import_array();
%}

// Interfaces
%include "../utils/PyreComponent.i"
%include "../problems/Physics.i"

%include "Material.i"
%include "Elasticity.i"
%include "RheologyElasticity.i"
%include "IsotropicLinearElasticity.i"
%include "IsotropicLinearMaxwell.i"
%include "IsotropicLinearGenMaxwell.i"
%include "IsotropicPowerLaw.i"
%include "IncompressibleElasticity.i"
%include "RheologyIncompressibleElasticity.i"
%include "IsotropicLinearIncompElasticity.i"
%include "Poroelasticity.i"
%include "RheologyPoroelasticity.i"
%include "IsotropicLinearPoroelasticity.i"

// End of file
