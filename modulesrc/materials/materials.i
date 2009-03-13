// -*- C++ -*-
//
// ======================================================================
//
//                           Brad T. Aagaard
//                        U.S. Geological Survey
//
// {LicenseText}
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
#include "pylith/materials/ElasticStrain1D.hh"
#include "pylith/materials/ElasticStress1D.hh"
#include "pylith/materials/ElasticPlaneStrain.hh"
#include "pylith/materials/ElasticPlaneStress.hh"
#include "pylith/materials/ElasticIsotropic3D.hh"

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
%include "../include/doublearray.i"

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
%include "ElasticStrain1D.i"
%include "ElasticStress1D.i"
%include "ElasticPlaneStrain.i"
%include "ElasticPlaneStress.i"
%include "ElasticIsotropic3D.i"


// End of file

