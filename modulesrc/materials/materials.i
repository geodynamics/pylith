// -*- C++ -*-
//
// ======================================================================
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2015 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

// SWIG interface
%module materials

// Header files for module C++ code
%{
#include "pylith/materials/materialsfwd.hh" // forward declarations

#include "spatialdata/spatialdb/spatialdbfwd.hh" // forward declarations
#include "spatialdata/units/unitsfwd.hh" // forward declarations

#include "pylith/materials/ElasticMaterial.hh"
#include "pylith/materials/ElasticPlaneStrain.hh"
#include "pylith/materials/ElasticPlaneStress.hh"
#include "pylith/materials/ElasticIsotropic3D.hh"
#include "pylith/materials/MaxwellIsotropic3D.hh"
#include "pylith/materials/MaxwellPlaneStrain.hh"
#include "pylith/materials/GenMaxwellIsotropic3D.hh"
#include "pylith/materials/GenMaxwellPlaneStrain.hh"
#include "pylith/materials/GenMaxwellQpQsIsotropic3D.hh"
#include "pylith/materials/PowerLaw3D.hh"
#include "pylith/materials/PowerLawPlaneStrain.hh"
#include "pylith/materials/DruckerPrager3D.hh"
#include "pylith/materials/DruckerPragerPlaneStrain.hh"

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
%include "Material.i"
%include "ElasticMaterial.i"
%include "ElasticPlaneStrain.i"
%include "ElasticPlaneStress.i"
%include "ElasticIsotropic3D.i"
%include "MaxwellIsotropic3D.i"
%include "MaxwellPlaneStrain.i"
%include "GenMaxwellIsotropic3D.i"
%include "GenMaxwellPlaneStrain.i"
%include "GenMaxwellQpQsIsotropic3D.i"
%include "PowerLaw3D.i"
%include "PowerLawPlaneStrain.i"
%include "DruckerPrager3D.i"
%include "DruckerPragerPlaneStrain.i"


// End of file

