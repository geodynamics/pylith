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
// Copyright (c) 2010-2012 University of California, Davis
//
// See COPYING for license information.
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
#include "pylith/topology/FieldsNew.hh"
#include "pylith/topology/SolutionFields.hh"
#include "pylith/topology/Jacobian.hh"
#include "pylith/topology/Distributor.hh"
#include "pylith/topology/RefineUniform.hh"
#include "pylith/topology/ReverseCuthillMcKee.hh"
%}

%include "exception.i"
%exception {
  try {
    $action
  } catch (const ALE::Exception& err) {
    SWIG_exception(SWIG_RuntimeError, err.message());
  } catch (const std::exception& err) {
    SWIG_exception(SWIG_RuntimeError, err.what());
  } // try/catch
 } // exception

%include "typemaps.i"
%include "../include/chararray.i"
%include "../include/submeshfield.i"
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
%include "Mesh.i"
%include "SubMesh.i"
%include "MeshOps.i"
%include "FieldBase.i"
%include "Field.i"
%include "Fields.i"
%include "FieldsNew.i"
%include "SolutionFields.i"
%include "Jacobian.i"
%include "Distributor.i"
%include "RefineUniform.i"
%include "ReverseCuthillMcKee.i"

// Template instatiation

%template(MeshField) pylith::topology::Field<pylith::topology::Mesh>;
%template(SubMeshField) pylith::topology::Field<pylith::topology::SubMesh>;
%template(MeshFields) pylith::topology::Fields<pylith::topology::Field<pylith::topology::Mesh> >;
%template(SubMeshFields) pylith::topology::Fields<pylith::topology::Field<pylith::topology::SubMesh> >;

%template(MeshFieldsNew) pylith::topology::FieldsNew<pylith::topology::Mesh>;
%template(SubMeshFieldsNew) pylith::topology::FieldsNew<pylith::topology::SubMesh>;

%extend pylith::topology::Field<pylith::topology::Mesh> {
  %template(createScatterMesh) createScatter<pylith::topology::Mesh>;
 }
%extend pylith::topology::Field<pylith::topology::Mesh> {
  %template(createScatterSubMesh) createScatter<pylith::topology::SubMesh>;
 }

// End of file

