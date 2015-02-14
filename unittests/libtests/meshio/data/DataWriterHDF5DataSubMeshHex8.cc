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

#include "DataWriterHDF5DataSubMeshHex8.hh"

#include <assert.h> // USES assert()

const char* pylith::meshio::DataWriterHDF5DataSubMeshHex8::_meshFilename = 
  "data/hex8.mesh";

const char* pylith::meshio::DataWriterHDF5DataSubMeshHex8::_bcLabel = 
  "top";

const char* pylith::meshio::DataWriterHDF5DataSubMeshHex8::_faultLabel = 0;
const int pylith::meshio::DataWriterHDF5DataSubMeshHex8::_faultId = 0;

const char* pylith::meshio::DataWriterHDF5DataSubMeshHex8::_timestepFilename = 
  "hex8_surf.h5";

const char* pylith::meshio::DataWriterHDF5DataSubMeshHex8::_vertexFilename = 
  "hex8_surf_vertex.h5";

const char* pylith::meshio::DataWriterHDF5DataSubMeshHex8::_cellFilename = 
  "hex8_surf_cell.h5";

const PylithScalar pylith::meshio::DataWriterHDF5DataSubMeshHex8::_time = 1.0;

const char* pylith::meshio::DataWriterHDF5DataSubMeshHex8::_timeFormat = 
  "%3.1f";

const int pylith::meshio::DataWriterHDF5DataSubMeshHex8::_numVertices = 12;

const pylith::meshio::DataWriterData::FieldStruct
pylith::meshio::DataWriterHDF5DataSubMeshHex8::_vertexFields[4] = {
  { "pressure", topology::FieldBase::SCALAR, 1 },
  { "displacement", topology::FieldBase::VECTOR, 3 },
  { "stress", topology::FieldBase::TENSOR, 6 },
  { "other", topology::FieldBase::OTHER, 2 },
};
const PylithScalar pylith::meshio::DataWriterHDF5DataSubMeshHex8::_vertexFieldScalar[12*1] = {
  2.1, 3.2, 4.3, 5.4, 6.5, 7.6, 8.7, 9.8, 10.9, 11.8, 12.7, 13.6
};
const PylithScalar pylith::meshio::DataWriterHDF5DataSubMeshHex8::_vertexFieldVector[12*3] = {
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
  7.9, 8.1, 9.2,
  10.3, 11.4, 12.5,
};
const PylithScalar pylith::meshio::DataWriterHDF5DataSubMeshHex8::_vertexFieldTensor[12*6] = {
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
};
const PylithScalar pylith::meshio::DataWriterHDF5DataSubMeshHex8::_vertexFieldOther[12*2] = {
  1.2, 2.3,
  3.4, 4.5,
  5.6, 6.7,
  7.8, 8.9,
  1.3, 2.4,
  3.5, 4.6,
  5.7, 6.8,
  7.9, 8.0,
  8.1, 8.2,
  9.2, 9.3,
  10.4, 10.5,
  11.5, 11.6,
};

const int pylith::meshio::DataWriterHDF5DataSubMeshHex8::_numCells = 2;

const pylith::meshio::DataWriterData::FieldStruct
pylith::meshio::DataWriterHDF5DataSubMeshHex8::_cellFields[4] = {
  { "pressure", topology::FieldBase::SCALAR, 1 },
  { "traction", topology::FieldBase::VECTOR, 3 },
  { "stress", topology::FieldBase::TENSOR, 6 },
  { "other", topology::FieldBase::OTHER, 3 },
};
const PylithScalar pylith::meshio::DataWriterHDF5DataSubMeshHex8::_cellFieldScalar[2*1] = {
  2.1, 3.2,
};
const PylithScalar pylith::meshio::DataWriterHDF5DataSubMeshHex8::_cellFieldVector[2*3] = {
  1.1, 2.2, 3.3,
  4.4, 5.5, 6.6,
};
const PylithScalar pylith::meshio::DataWriterHDF5DataSubMeshHex8::_cellFieldTensor[2*6] = {
  1.2, 2.3, 3.4, 4.5, 5.6, 6.7,
  7.8, 8.9, 9.0, 10.1, 11.2, 12.3,
};
const PylithScalar pylith::meshio::DataWriterHDF5DataSubMeshHex8::_cellFieldOther[2*3] = {
  1.2, 2.3, 3.4,
  7.8, 8.9, 9.0,
};

pylith::meshio::DataWriterHDF5DataSubMeshHex8::DataWriterHDF5DataSubMeshHex8(void)
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

pylith::meshio::DataWriterHDF5DataSubMeshHex8::~DataWriterHDF5DataSubMeshHex8(void)
{}


// End of file
