// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================
//

// Define the materialscontrib SWIG interface module.

// Set module name
%module materialscontrib

// Header files for module C++ code.
%{
#include "pylith/materials/materialsfwd.hh" // forward declarations

#include "spatialdata/spatialdb/spatialdbfwd.hh" // forward declarations
#include "spatialdata/units/unitsfwd.hh" // forward declarations

#include "PlaneStrainState.hh"

#include "pylith/utils/arrayfwd.hh"
    %
}

// Convert standard C++ exceptions to Python exceptions.
%include "exception.i"
%exception {
    try {
        $ action
    } catch (const std::exception& err) {
        SWIG_exception (SWIG_RuntimeError, err.what ());
    } // try/catch
} // exception

%include "typemaps.i"
%include "include/scalartypemaps.i"

// Numpy interface stuff
%{
#define SWIG_FILE_WITH_INIT
    %
}
%include "include/numpy.i"
%init %{
    import_array();
    %
}

// Interface files.
%include "materials/Material.i"
%include "materials/ElasticMaterial.i"
%include "PlaneStrainState.i"

// End of file
