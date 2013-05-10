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
// Copyright (c) 2010-2013 University of California, Davis
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
#if defined(ENABLE_CUBIT)
#include "pylith/meshio/MeshIOCubit.hh"
#endif

#include "pylith/meshio/VertexFilter.hh"
#include "pylith/meshio/VertexFilterVecNorm.hh"
#include "pylith/meshio/CellFilter.hh"
#include "pylith/meshio/CellFilterAvg.hh"
#include "pylith/meshio/DataWriter.hh"
#include "pylith/meshio/DataWriterVTK.hh"
#include "pylith/meshio/OutputManager.hh"
#include "pylith/meshio/OutputSolnSubset.hh"
#include "pylith/meshio/OutputSolnPoints.hh"
#if defined(ENABLE_HDF5)
#include "pylith/meshio/DataWriterHDF5.hh"
#include "pylith/meshio/DataWriterHDF5Ext.hh"
#include "pylith/meshio/Xdmf.hh"
#endif

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
%include "MeshIOObj.i"
%include "MeshIOAscii.i"
%include "MeshIOLagrit.i"
#if defined(ENABLE_CUBIT)
%include "MeshIOCubit.i"
#endif

%include "VertexFilter.i"
%include "VertexFilterVecNorm.i"
%include "CellFilter.i"
%include "CellFilterAvg.i"
%include "DataWriter.i"
%include "DataWriterVTK.i"
%include "OutputManager.i"
%include "OutputSolnSubset.i"
%include "OutputSolnPoints.i"
#if defined(ENABLE_HDF5)
%include "DataWriterHDF5.i"
%include "DataWriterHDF5Ext.i"
%include "Xdmf.i"
#endif

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
%template(PointsDataWriterVTK) pylith::meshio::DataWriterVTK<pylith::topology::Mesh, pylith::topology::Field<pylith::topology::Mesh> >;

#if defined(ENABLE_HDF5)
%template(MeshDataWriterHDF5) pylith::meshio::DataWriterHDF5<pylith::topology::Mesh, pylith::topology::Field<pylith::topology::Mesh> >;
%template(SubMeshDataWriterHDF5) pylith::meshio::DataWriterHDF5<pylith::topology::SubMesh, pylith::topology::Field<pylith::topology::Mesh> >;
%template(SubSubMeshDataWriterHDF5) pylith::meshio::DataWriterHDF5<pylith::topology::SubMesh, pylith::topology::Field<pylith::topology::SubMesh> >;

%template(MeshDataWriterHDF5Ext) pylith::meshio::DataWriterHDF5Ext<pylith::topology::Mesh, pylith::topology::Field<pylith::topology::Mesh> >;
%template(SubMeshDataWriterHDF5Ext) pylith::meshio::DataWriterHDF5Ext<pylith::topology::SubMesh, pylith::topology::Field<pylith::topology::Mesh> >;
%template(SubSubMeshDataWriterHDF5Ext) pylith::meshio::DataWriterHDF5Ext<pylith::topology::SubMesh, pylith::topology::Field<pylith::topology::SubMesh> >;
#endif

%template(MeshOutputManager) pylith::meshio::OutputManager<pylith::topology::Mesh, pylith::topology::Field<pylith::topology::Mesh> >;
%template(SubMeshOutputManager) pylith::meshio::OutputManager<pylith::topology::SubMesh, pylith::topology::Field<pylith::topology::SubMesh> >;

// End of file
