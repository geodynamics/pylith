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

#include "DataWriterVTKDataBCMeshHex8.hh"

#include <assert.h> // USES assert()

const char* pylith::meshio::DataWriterVTKDataBCMeshHex8::_meshFilename = 
  "data/hex8.mesh";

const char* pylith::meshio::DataWriterVTKDataBCMeshHex8::_bcLabel = 
  "top";

const char* pylith::meshio::DataWriterVTKDataBCMeshHex8::_faultLabel = 
  "fault";
const int pylith::meshio::DataWriterVTKDataBCMeshHex8::_faultId = 100;

const char* pylith::meshio::DataWriterVTKDataBCMeshHex8::_timestepFilename = 
  "hex8_bc.vtk";

const char* pylith::meshio::DataWriterVTKDataBCMeshHex8::_vertexFilename = 
  "hex8_bc_vertex.vtk";

const char* pylith::meshio::DataWriterVTKDataBCMeshHex8::_cellFilename = 
  "hex8_bc_cell.vtk";

const PylithScalar pylith::meshio::DataWriterVTKDataBCMeshHex8::_time = 1.0;

const char* pylith::meshio::DataWriterVTKDataBCMeshHex8::_timeFormat = 
  "%3.1f";

const int pylith::meshio::DataWriterVTKDataBCMeshHex8::_numVertices = 8;

const pylith::meshio::DataWriterData::FieldStruct
pylith::meshio::DataWriterVTKDataBCMeshHex8::_vertexFields[4] = {
  { "pressure", topology::FieldBase::SCALAR, 1 },
  { "displacement", topology::FieldBase::VECTOR, 3 },
  { "stress", topology::FieldBase::TENSOR, 6 },
  { "other", topology::FieldBase::OTHER, 2 },
};
const PylithScalar pylith::meshio::DataWriterVTKDataBCMeshHex8::_vertexFieldScalar[8*1] = {
  2.1, 3.2, 4.3, 5.4, 6.5, 7.6, 8.7, 9.8,
};
const PylithScalar pylith::meshio::DataWriterVTKDataBCMeshHex8::_vertexFieldVector[8*3] = {
  1.1, 2.2, 3.3,
  4.4, 5.5, 6.6,
  7.7, 8.8, 9.9,
  10.1, 11.2, 12.3,
  1.2, 2.3, 3.4,
  4.5, 5.6, 6.7,
  7.8, 8.9, 9.0,
  10.2, 11.3, 12.4,
};
const PylithScalar pylith::meshio::DataWriterVTKDataBCMeshHex8::_vertexFieldTensor[8*6] = {
  1.1, 1.2, 1.3, 1.4, 1.5, 1.6,
  2.1, 2.2, 2.3, 2.4, 2.5, 2.6,
  3.1, 3.2, 3.3, 3.4, 3.5, 3.6,
  4.1, 4.2, 4.3, 4.4, 4.5, 4.6,
  5.1, 5.2, 5.3, 5.4, 5.5, 5.6,
  6.1, 6.2, 6.3, 6.4, 6.5, 6.6,
  7.1, 7.2, 7.3, 7.4, 7.5, 7.6,
  8.1, 8.2, 8.3, 8.4, 8.5, 8.6,
};
const PylithScalar pylith::meshio::DataWriterVTKDataBCMeshHex8::_vertexFieldOther[8*2] = {
  1.2, 2.3,
  3.4, 4.5,
  5.6, 6.7,
  7.8, 8.9,
  1.3, 2.4,
  3.5, 4.6,
  5.7, 6.8,
  7.9, 8.0,
};

const int pylith::meshio::DataWriterVTKDataBCMeshHex8::_numCells = 2;

const pylith::meshio::DataWriterData::FieldStruct
pylith::meshio::DataWriterVTKDataBCMeshHex8::_cellFields[4] = {
  { "pressure", topology::FieldBase::SCALAR, 1 },
  { "traction", topology::FieldBase::VECTOR, 3 },
  { "stress", topology::FieldBase::TENSOR, 6 },
  { "other", topology::FieldBase::OTHER, 2 },
};
const PylithScalar pylith::meshio::DataWriterVTKDataBCMeshHex8::_cellFieldScalar[2*1] = {
  2.1, 3.2,
};
const PylithScalar pylith::meshio::DataWriterVTKDataBCMeshHex8::_cellFieldVector[2*3] = {
  1.1, 2.2, 3.3,
  4.4, 5.5, 6.6,
};
const PylithScalar pylith::meshio::DataWriterVTKDataBCMeshHex8::_cellFieldTensor[2*6] = {
  1.2, 2.3, 3.4, 4.5, 5.6, 6.7,
  7.8, 8.9, 9.0, 10.1, 11.2, 12.3,
};
const PylithScalar pylith::meshio::DataWriterVTKDataBCMeshHex8::_cellFieldOther[2*2] = {
  1.1, 2.2,
  5.5, 6.6,
};

pylith::meshio::DataWriterVTKDataBCMeshHex8::DataWriterVTKDataBCMeshHex8(void)
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

pylith::meshio::DataWriterVTKDataBCMeshHex8::~DataWriterVTKDataBCMeshHex8(void)
{}


// End of file
