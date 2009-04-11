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
#include "pylith/topology/MeshOps.hh"
#include "pylith/topology/FieldBase.hh"
#include "pylith/topology/Field.hh"
#include "pylith/topology/Fields.hh"
#include "pylith/topology/SolutionFields.hh"
#include "pylith/topology/Jacobian.hh"
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
%include "../include/chararray.i"
%include "../include/submeshfield.i"

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
%include "MeshOps.i"
%include "FieldBase.i"
%include "Field.i"
%include "Fields.i"
%include "SolutionFields.i"
%include "Jacobian.i"

// Template instatiation
%template(MeshField) pylith::topology::Field<pylith::topology::Mesh>;
%template(SubMeshField) pylith::topology::Field<pylith::topology::SubMesh>;
%template(MeshFields) pylith::topology::Fields<pylith::topology::Field<pylith::topology::Mesh> >;
%template(SubMeshFields) pylith::topology::Fields<pylith::topology::Field<pylith::topology::SubMesh> >;

// End of file

