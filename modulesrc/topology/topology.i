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
%module topology

// Header files for module C++ code
%{
#include "pylith/topology/Mesh.hh"
#include "pylith/topology/SubMesh.hh"
// #include "pylith/topology/MeshOps.hh"
// #include "pylith/topology/FieldBase.hh"
#include "pylith/topology/Field.hh"
// #include "pylith/topology/FieldSubMesh.hh"
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

// Numpy interface stuff
%{
#define SWIG_FILE_WITH_INIT
%}
%include "../include/numpy.i"
%init %{
import_array();
%}

// Interfaces
%include "Mesh.i"
%include "SubMesh.i"
// %include "MeshOps.i"
// %include "FieldBase.i"
%include "Field.i"
// %include "FieldSubMesh.i"


// Template instatiation
%template(MeshField) pylith::topology::Field<pylith::topology::Mesh>;
%template(SubMeshField) pylith::topology::Field<pylith::topology::SubMesh>;

// End of file

