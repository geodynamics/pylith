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

#include "DataWriterVTKDataMatMeshHex8.hh"

#include <assert.h> // USES assert()

const char* pylith::meshio::DataWriterVTKDataMatMeshHex8::_meshFilename = 
  "data/hex8.mesh";

const char* pylith::meshio::DataWriterVTKDataMatMeshHex8::_cellsLabel = 
  "material-id";
const int pylith::meshio::DataWriterVTKDataMatMeshHex8::_labelId = 0;

const char* pylith::meshio::DataWriterVTKDataMatMeshHex8::_faultLabel = 
  "fault";
const int pylith::meshio::DataWriterVTKDataMatMeshHex8::_faultId = 100;

const char* pylith::meshio::DataWriterVTKDataMatMeshHex8::_timestepFilename = 
  "hex8_mat.vtk";

const char* pylith::meshio::DataWriterVTKDataMatMeshHex8::_vertexFilename = 
  "hex8_mat_vertex.vtk";

const char* pylith::meshio::DataWriterVTKDataMatMeshHex8::_cellFilename = 
  "hex8_mat_cell.vtk";

const PylithScalar pylith::meshio::DataWriterVTKDataMatMeshHex8::_time = 1.0;

const char* pylith::meshio::DataWriterVTKDataMatMeshHex8::_timeFormat = 
  "%3.1f";

const int pylith::meshio::DataWriterVTKDataMatMeshHex8::_numVertices = 16;

const pylith::meshio::DataWriterData::FieldStruct
pylith::meshio::DataWriterVTKDataMatMeshHex8::_vertexFields[4] = {
  { "pressure", topology::FieldBase::SCALAR, 1 },
  { "displacement", topology::FieldBase::VECTOR, 3 },
  { "stress", topology::FieldBase::TENSOR, 6 },
  { "other", topology::FieldBase::OTHER, 2 },
};
const PylithScalar pylith::meshio::DataWriterVTKDataMatMeshHex8::_vertexFieldScalar[16*1] = {
  2.1, 3.2, 4.3, 5.4, 6.5, 7.6, 8.7, 9.8,
  10.0, 12.1, 11.1, 13.1, 14.1, 15.1, 16.1, 17.1,
};
const PylithScalar pylith::meshio::DataWriterVTKDataMatMeshHex8::_vertexFieldVector[16*3] = {
  1.1, 2.2, 3.3,
  4.4, 5.5, 6.6,
  7.7, 8.8, 9.9,
  10.1, 11.2, 12.3,
  1.2, 2.3, 3.4,
  4.5, 5.6, 6.7,
  7.8, 8.9, 9.0,
  10.2, 11.3, 12.4,
  1.3, 2.4, 3.5,
  4.6, 5.7, 6.8,
  7.9, 8.0, 9.1,
  10.2, 11.3, 12.4,
  13.5, 14.6, 15.7,
  16.8, 17.9, 18.1,
  19.2, 20.3, 21.4,
  22.5, 23.6, 24.7,
};
const PylithScalar pylith::meshio::DataWriterVTKDataMatMeshHex8::_vertexFieldTensor[16*6] = {
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
  12.1, 12.2, 12.3, 12.4, 12.5, 12.6,
  13.1, 13.2, 13.3, 13.4, 13.5, 13.6,
  14.1, 14.2, 14.3, 14.4, 14.5, 14.6,
  15.1, 15.2, 15.3, 15.4, 15.5, 15.6,
  16.1, 16.2, 16.3, 16.4, 16.5, 16.6,
};
const PylithScalar pylith::meshio::DataWriterVTKDataMatMeshHex8::_vertexFieldOther[16*2] = {
  1.2, 2.3,
  3.4, 4.5,
  5.6, 6.7,
  7.8, 8.9,
  1.3, 2.4,
  3.5, 4.6,
  5.7, 6.8,
  7.9, 8.0,
  1.3, 2.4,
  3.5, 4.6,
  5.7, 6.8,
  8.0, 1.4,
  2.5, 3.6,
  4.8, 1.5,
  2.6, 3.7,
  4.8, 5.9,
};

const int pylith::meshio::DataWriterVTKDataMatMeshHex8::_numCells = 1;

const pylith::meshio::DataWriterData::FieldStruct
pylith::meshio::DataWriterVTKDataMatMeshHex8::_cellFields[4] = {
  { "pressure", topology::FieldBase::SCALAR, 1 },
  { "traction", topology::FieldBase::VECTOR, 3 },
  { "stress", topology::FieldBase::TENSOR, 6 },
  { "other", topology::FieldBase::OTHER, 2 },
};
const PylithScalar pylith::meshio::DataWriterVTKDataMatMeshHex8::_cellFieldScalar[1*1] = {
  2.1,
};
const PylithScalar pylith::meshio::DataWriterVTKDataMatMeshHex8::_cellFieldVector[1*3] = {
  1.1, 2.2, 3.3,
};
const PylithScalar pylith::meshio::DataWriterVTKDataMatMeshHex8::_cellFieldTensor[1*6] = {
  1.2, 2.3, 3.4, 4.5, 5.6, 6.7,
};
const PylithScalar pylith::meshio::DataWriterVTKDataMatMeshHex8::_cellFieldOther[1*2] = {
  1.2, 2.3,
};

pylith::meshio::DataWriterVTKDataMatMeshHex8::DataWriterVTKDataMatMeshHex8(void)
{ // constructor
  meshFilename = const_cast<char*>(_meshFilename);
  cellsLabel = const_cast<char*>(_cellsLabel);
  labelId = _labelId;
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

pylith::meshio::DataWriterVTKDataMatMeshHex8::~DataWriterVTKDataMatMeshHex8(void)
{}


// End of file
