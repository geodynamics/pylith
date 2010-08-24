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
// Copyright (c) 2010 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

// SWIG interface
%module meshio

// Header files for module C++ code
%{
#include "pylith/meshio/MeshIO.hh"
#include "pylith/meshio/MeshIOAscii.hh"
#include "pylith/meshio/MeshIOLagrit.hh"
#include "pylith/meshio/MeshIOSieve.hh"
#ifdef ENABLE_CUBIT
#include "pylith/meshio/MeshIOCubit.hh"
#endif

#include "pylith/meshio/VertexFilter.hh"
#include "pylith/meshio/VertexFilterVecNorm.hh"
#include "pylith/meshio/CellFilter.hh"
#include "pylith/meshio/CellFilterAvg.hh"
#include "pylith/meshio/DataWriter.hh"
#include "pylith/meshio/DataWriterVTK.hh"
#include "pylith/meshio/DataWriterHDF5.hh"
#include "pylith/meshio/OutputManager.hh"
#include "pylith/meshio/OutputSolnSubset.hh"

#include "pylith/utils/arrayfwd.hh"
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

// Interfaces
%include "MeshIOObj.i"
%include "MeshIOAscii.i"
%include "MeshIOLagrit.i"
%include "MeshIOSieve.i"
#ifdef ENABLE_CUBIT
%include "MeshIOCubit.i"
#endif

%include "VertexFilter.i"
%include "VertexFilterVecNorm.i"
%include "CellFilter.i"
%include "CellFilterAvg.i"
%include "DataWriter.i"
%include "DataWriterVTK.i"
%include "DataWriterHDF5.i"
%include "OutputManager.i"
%include "OutputSolnSubset.i"

// Template instatiation
%template(MeshVertexFilter) pylith::meshio::VertexFilter<pylith::topology::Field<pylith::topology::Mesh> >;
%template(SubMeshVertexFilter) pylith::meshio::VertexFilter<pylith::topology::Field<pylith::topology::SubMesh> >;
%template(MeshVertexFilterVecNorm) pylith::meshio::VertexFilterVecNorm<pylith::topology::Field<pylith::topology::Mesh> >;
%template(SubMeshVertexFilterVecNorm) pylith::meshio::VertexFilterVecNorm<pylith::topology::Field<pylith::topology::SubMesh> >;

%template(MeshCellFilter) pylith::meshio::CellFilter<pylith::topology::Mesh, pylith::topology::Field<pylith::topology::Mesh> >;
%template(SubMeshCellFilter) pylith::meshio::CellFilter<pylith::topology::SubMesh, pylith::topology::Field<pylith::topology::SubMesh> >;
%template(MeshCellFilterAvg) pylith::meshio::CellFilterAvg<pylith::topology::Mesh, pylith::topology::Field<pylith::topology::Mesh> >;
%template(SubMeshCellFilterAvg) pylith::meshio::CellFilterAvg<pylith::topology::SubMesh, pylith::topology::Field<pylith::topology::SubMesh> >;

%template(MeshDataWriter) pylith::meshio::DataWriter<pylith::topology::Mesh, pylith::topology::Field<pylith::topology::Mesh> >;
%template(SubMeshDataWriter) pylith::meshio::DataWriter<pylith::topology::SubMesh, pylith::topology::Field<pylith::topology::Mesh> >;
%template(SubSubMeshDataWriter) pylith::meshio::DataWriter<pylith::topology::SubMesh, pylith::topology::Field<pylith::topology::SubMesh> >;

%template(MeshDataWriterVTK) pylith::meshio::DataWriterVTK<pylith::topology::Mesh, pylith::topology::Field<pylith::topology::Mesh> >;
%template(SubMeshDataWriterVTK) pylith::meshio::DataWriterVTK<pylith::topology::SubMesh, pylith::topology::Field<pylith::topology::Mesh> >;
%template(SubSubMeshDataWriterVTK) pylith::meshio::DataWriterVTK<pylith::topology::SubMesh, pylith::topology::Field<pylith::topology::SubMesh> >;

%template(MeshDataWriterHDF5) pylith::meshio::DataWriterHDF5<pylith::topology::Mesh, pylith::topology::Field<pylith::topology::Mesh> >;
%template(SubMeshDataWriterHDF5) pylith::meshio::DataWriterHDF5<pylith::topology::SubMesh, pylith::topology::Field<pylith::topology::Mesh> >;
%template(SubSubMeshDataWriterHDF5) pylith::meshio::DataWriterHDF5<pylith::topology::SubMesh, pylith::topology::Field<pylith::topology::SubMesh> >;

%template(MeshOutputManager) pylith::meshio::OutputManager<pylith::topology::Mesh, pylith::topology::Field<pylith::topology::Mesh> >;
%template(SubMeshOutputManager) pylith::meshio::OutputManager<pylith::topology::SubMesh, pylith::topology::Field<pylith::topology::SubMesh> >;

// End of file

