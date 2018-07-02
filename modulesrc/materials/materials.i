// ======================================================================
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2017 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

// SWIG interface
%module materials

// Header files for module C++ code
%{
#include "pylith/materials/Material.hh"
#include "pylith/materials/IsotropicLinearElasticityPlaneStrain.hh"
#include "pylith/materials/IsotropicLinearIncompElasticityPlaneStrain.hh"
#include "pylith/materials/IsotropicLinearMaxwellPlaneStrain.hh"
#include "pylith/materials/IsotropicLinearGenMaxwellPlaneStrain.hh"
#include "pylith/materials/IsotropicLinearElasticity3D.hh"
#include "pylith/materials/IsotropicLinearMaxwell3D.hh"
#include "pylith/materials/IsotropicLinearGenMaxwell3D.hh"

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
%include "../feassemble/ObservedComponent.i"
%include "../feassemble/IntegratorPointwise.i"

%include "Material.i"
%include "IsotropicLinearElasticityPlaneStrain.i"
%include "IsotropicLinearIncompElasticityPlaneStrain.i"
%include "IsotropicLinearMaxwellPlaneStrain.i"
%include "IsotropicLinearGenMaxwellPlaneStrain.i"
%include "IsotropicLinearElasticity3D.i"
%include "IsotropicLinearMaxwell3D.i"
%include "IsotropicLinearGenMaxwell3D.i"


// End of file
