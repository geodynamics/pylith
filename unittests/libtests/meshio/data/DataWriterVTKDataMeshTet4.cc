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

#include "DataWriterVTKDataMeshTet4.hh"

#include <assert.h> // USES assert()

const char* pylith::meshio::DataWriterVTKDataMeshTet4::_meshFilename = 
  "data/tet4.mesh";

const char* pylith::meshio::DataWriterVTKDataMeshTet4::_faultLabel = 
  "fault";
const int pylith::meshio::DataWriterVTKDataMeshTet4::_faultId = 100;

const char* pylith::meshio::DataWriterVTKDataMeshTet4::_timestepFilename = 
  "tet4.vtk";

const char* pylith::meshio::DataWriterVTKDataMeshTet4::_vertexFilename = 
  "tet4_vertex.vtk";

const char* pylith::meshio::DataWriterVTKDataMeshTet4::_cellFilename = 
  "tet4_cell.vtk";

const PylithScalar pylith::meshio::DataWriterVTKDataMeshTet4::_time = 1.0;

const char* pylith::meshio::DataWriterVTKDataMeshTet4::_timeFormat = 
  "%3.1f";

const int pylith::meshio::DataWriterVTKDataMeshTet4::_numVertices = 11;

const pylith::meshio::DataWriterData::FieldStruct
pylith::meshio::DataWriterVTKDataMeshTet4::_vertexFields[4] = {
  { "pressure", topology::FieldBase::SCALAR, 1 },
  { "displacement", topology::FieldBase::VECTOR, 3 },
  { "stress", topology::FieldBase::TENSOR, 6 },
  { "other", topology::FieldBase::OTHER, 2 },
};
const PylithScalar pylith::meshio::DataWriterVTKDataMeshTet4::_vertexFieldScalar[11*1] = {
  2.1, 3.2, 4.3, 5.4, 6.5, 7.6, 8.7, 9.8, 10.9, 11.0, 12.1,
};
const PylithScalar pylith::meshio::DataWriterVTKDataMeshTet4::_vertexFieldVector[11*3] = {
  1.1, 2.2, 3.3,
  4.4, 5.5, 6.6,
  7.7, 8.8, 9.9,
  10.0, 11.1, 12.2,
  13.3, 14.4, 15.5,
  16.6, 17.7, 18.8,
  19.9, 20.0, 21.1,
  22.2, 23.3, 24.4,
  25.5, 26.6, 27.7,
  28.8, 29.9, 30.0,
  31.1, 32.2, 33.3,
};
const PylithScalar pylith::meshio::DataWriterVTKDataMeshTet4::_vertexFieldTensor[11*6] = {
  1.1, 1.2, 1.3, 1.4, 1.5, 1.6,
  2.1, 2.2, 2.3, 2.4, 2.5, 2.6,
  3.1, 3.2, 3.3, 3.4, 3.5, 3.6,
  4.1, 4.2, 4.3, 4.4, 4.5, 4.6,
  5.1, 5.2, 5.3, 5.4, 5.5, 5.6,
  6.1, 6.2, 6.3, 6.4, 6.5, 6.6,
  7.1, 7.2, 7.3, 7.4, 7.5, 7.6,
  8.1, 8.2, 8.3, 8.4, 8.5, 8.6,
  9.1, 9.2, 9.3, 9.4, 9.5, 9.6,
  10.1, 10.2, 10.3, 10.4, 10.5, 10.6,
  11.1, 11.2, 11.3, 11.4, 11.5, 11.6,
};
const PylithScalar pylith::meshio::DataWriterVTKDataMeshTet4::_vertexFieldOther[11*2] = {
  1.2, 2.3,
  3.4, 4.5,
  5.6, 6.7,
  7.8, 8.9,
  9.0, 10.1,
  11.2, 12.3,
  13.4, 14.5,
  15.6, 16.7,
  17.8, 18.9,
  19.0, 20.1,
  21.2, 22.3,
};

const int pylith::meshio::DataWriterVTKDataMeshTet4::_numCells = 3;

const pylith::meshio::DataWriterData::FieldStruct
pylith::meshio::DataWriterVTKDataMeshTet4::_cellFields[4] = {
  { "pressure", topology::FieldBase::SCALAR, 1 },
  { "traction", topology::FieldBase::VECTOR, 3 },
  { "stress", topology::FieldBase::TENSOR, 6 },
  { "other", topology::FieldBase::OTHER, 4 },
};
const PylithScalar pylith::meshio::DataWriterVTKDataMeshTet4::_cellFieldScalar[3*1] = {
  2.1, 3.2, 4.3
};
const PylithScalar pylith::meshio::DataWriterVTKDataMeshTet4::_cellFieldVector[3*3] = {
  1.1, 2.2, 3.3,
  4.4, 5.5, 6.6,
  7.7, 8.8, 9.9
};
const PylithScalar pylith::meshio::DataWriterVTKDataMeshTet4::_cellFieldTensor[3*6] = {
  1.2, 2.3, 3.4, 4.5, 5.6, 6.7,
  7.8, 8.9, 9.0, 10.1, 11.2, 12.3,
  13.4, 14.5, 15.6, 16.7, 17.8, 18.9
};
const PylithScalar pylith::meshio::DataWriterVTKDataMeshTet4::_cellFieldOther[3*4] = {
  1.2, 2.3, 3.4, 4.5,
  7.8, 8.9, 9.0, 10.1,
  13.4, 14.5, 15.6, 16.7,
};

pylith::meshio::DataWriterVTKDataMeshTet4::DataWriterVTKDataMeshTet4(void)
{ // constructor
  meshFilename = const_cast<char*>(_meshFilename);
  faultLabel = const_cast<char*>(_faultLabel);
  faultId = _faultId;

  timestepFilename = const_cast<char*>(_timestepFilename);
  vertexFilename = const_cast<char*>(_vertexFilename);
  cellFilename = const_cast<char*>(_cellFilename);

  time = _time;
  timeFormat = const_cast<char*>(_timeFormat);
  
  numVertices = _numVertices;
  assert(DataWriterData::numVertexFields == numVertexFields);
  vertexFieldsInfo = const_cast<DataWriterData::FieldStruct*>(_vertexFields);
  vertexFields[0] = const_cast<PylithScalar*>(_vertexFieldScalar);
  vertexFields[1] = const_cast<PylithScalar*>(_vertexFieldVector);
  vertexFields[2] = const_cast<PylithScalar*>(_vertexFieldTensor);
  vertexFields[3] = const_cast<PylithScalar*>(_vertexFieldOther);

  numCells = _numCells;
  assert(DataWriterData::numCellFields == numCellFields);
  cellFieldsInfo = const_cast<DataWriterData::FieldStruct*>(_cellFields);
  cellFields[0] = const_cast<PylithScalar*>(_cellFieldScalar);
  cellFields[1] = const_cast<PylithScalar*>(_cellFieldVector);
  cellFields[2] = const_cast<PylithScalar*>(_cellFieldTensor);
  cellFields[3] = const_cast<PylithScalar*>(_cellFieldOther);
} // constructor

pylith::meshio::DataWriterVTKDataMeshTet4::~DataWriterVTKDataMeshTet4(void)
{}


// End of file
