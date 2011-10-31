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
// Copyright (c) 2010-2011 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

#include "DataWriterVTKDataSubMeshTet4.hh"

#include <assert.h> // USES assert()

const char* pylith::meshio::DataWriterVTKDataSubMeshTet4::_meshFilename = 
  "data/tet4.mesh";

const char* pylith::meshio::DataWriterVTKDataSubMeshTet4::_bcLabel = 
  "boundary";

const char* pylith::meshio::DataWriterVTKDataSubMeshTet4::_faultLabel =  0;
const int pylith::meshio::DataWriterVTKDataSubMeshTet4::_faultId = 0;

const char* pylith::meshio::DataWriterVTKDataSubMeshTet4::_timestepFilename = 
  "tet4_surf.vtk";

const char* pylith::meshio::DataWriterVTKDataSubMeshTet4::_vertexFilename = 
  "tet4_surf_vertex.vtk";

const char* pylith::meshio::DataWriterVTKDataSubMeshTet4::_cellFilename = 
  "tet4_surf_cell.vtk";

const PylithScalar pylith::meshio::DataWriterVTKDataSubMeshTet4::_time = 1.0;

const char* pylith::meshio::DataWriterVTKDataSubMeshTet4::_timeFormat = 
  "%3.1f";

const int pylith::meshio::DataWriterVTKDataSubMeshTet4::_numVertexFields = 3;
const int pylith::meshio::DataWriterVTKDataSubMeshTet4::_numVertices = 5;

const pylith::meshio::DataWriterData::FieldStruct
pylith::meshio::DataWriterVTKDataSubMeshTet4::_vertexFields[] = {
  { "displacements", topology::FieldBase::VECTOR, 3 },
  { "pressure", topology::FieldBase::SCALAR, 1 },
  { "other", topology::FieldBase::OTHER, 2 },
};
const PylithScalar pylith::meshio::DataWriterVTKDataSubMeshTet4::_vertexField0[] = {
  1.1, 2.2, 3.3,
  4.4, 5.5, 6.6,
  7.7, 8.8, 9.9,
  10.0, 11.1, 12.2,
  13.3, 14.4, 15.5,
};
const PylithScalar pylith::meshio::DataWriterVTKDataSubMeshTet4::_vertexField1[] = {
  2.1, 3.2, 4.3, 5.4, 6.5,
};
const PylithScalar pylith::meshio::DataWriterVTKDataSubMeshTet4::_vertexField2[] = {
  1.2, 2.3,
  3.4, 4.5,
  5.6, 6.7,
  7.8, 8.9,
  9.0, 10.1,
};

const int pylith::meshio::DataWriterVTKDataSubMeshTet4::_numCellFields = 3;
const int pylith::meshio::DataWriterVTKDataSubMeshTet4::_numCells = 2;

const pylith::meshio::DataWriterData::FieldStruct
pylith::meshio::DataWriterVTKDataSubMeshTet4::_cellFields[] = {
  { "traction", topology::FieldBase::VECTOR, 3 },
  { "pressure", topology::FieldBase::SCALAR, 1 },
  { "other", topology::FieldBase::TENSOR, 6 },
};
const PylithScalar pylith::meshio::DataWriterVTKDataSubMeshTet4::_cellField0[] = {
  1.1, 2.2, 3.3,
  4.4, 5.5, 6.6,
};
const PylithScalar pylith::meshio::DataWriterVTKDataSubMeshTet4::_cellField1[] = {
  2.1, 3.2
};
const PylithScalar pylith::meshio::DataWriterVTKDataSubMeshTet4::_cellField2[] = {
  1.2, 2.3, 3.4, 4.5, 5.6, 6.7,
  7.8, 8.9, 9.0, 10.1, 11.2, 12.3,
};

pylith::meshio::DataWriterVTKDataSubMeshTet4::DataWriterVTKDataSubMeshTet4(void)
{ // constructor
  meshFilename = const_cast<char*>(_meshFilename);
  bcLabel = const_cast<char*>(_bcLabel);
  faultLabel = const_cast<char*>(_faultLabel);
  faultId = _faultId;

  timestepFilename = const_cast<char*>(_timestepFilename);
  vertexFilename = const_cast<char*>(_vertexFilename);
  cellFilename = const_cast<char*>(_cellFilename);

  time = _time;
  timeFormat = const_cast<char*>(_timeFormat);
  
  numVertexFields = _numVertexFields;
  assert(3 == numVertexFields);
  numVertices = _numVertices;
  vertexFieldsInfo = const_cast<DataWriterData::FieldStruct*>(_vertexFields);
  vertexFields[0] = const_cast<PylithScalar*>(_vertexField0);
  vertexFields[1] = const_cast<PylithScalar*>(_vertexField1);
  vertexFields[2] = const_cast<PylithScalar*>(_vertexField2);

  numCellFields = _numCellFields;
  assert(3 == numCellFields);
  numCells = _numCells;
  cellFieldsInfo = const_cast<DataWriterData::FieldStruct*>(_cellFields);
  cellFields[0] = const_cast<PylithScalar*>(_cellField0);
  cellFields[1] = const_cast<PylithScalar*>(_cellField1);
  cellFields[2] = const_cast<PylithScalar*>(_cellField2);
} // constructor

pylith::meshio::DataWriterVTKDataSubMeshTet4::~DataWriterVTKDataSubMeshTet4(void)
{}


// End of file
