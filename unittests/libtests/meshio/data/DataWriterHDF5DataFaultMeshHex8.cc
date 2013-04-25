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

#include "DataWriterHDF5DataFaultMeshHex8.hh"

#include <assert.h> // USES assert()

const char* pylith::meshio::DataWriterHDF5DataFaultMeshHex8::_meshFilename = 
  "data/hex8.mesh";

const char* pylith::meshio::DataWriterHDF5DataFaultMeshHex8::_faultLabel = 
  "fault";
const int pylith::meshio::DataWriterHDF5DataFaultMeshHex8::_faultId = 100;

const char* pylith::meshio::DataWriterHDF5DataFaultMeshHex8::_timestepFilename = 
  "hex8_fault.h5";

const char* pylith::meshio::DataWriterHDF5DataFaultMeshHex8::_vertexFilename = 
  "hex8_fault_vertex.h5";

const char* pylith::meshio::DataWriterHDF5DataFaultMeshHex8::_cellFilename = 
  "hex8_fault_cell.h5";

const PylithScalar pylith::meshio::DataWriterHDF5DataFaultMeshHex8::_time = 1.0;

const char* pylith::meshio::DataWriterHDF5DataFaultMeshHex8::_timeFormat = 
  "%3.1f";

const int pylith::meshio::DataWriterHDF5DataFaultMeshHex8::_numVertices = 4;

const pylith::meshio::DataWriterData::FieldStruct
pylith::meshio::DataWriterHDF5DataFaultMeshHex8::_vertexFields[4] = {
  { "pressure", topology::FieldBase::SCALAR, 1 },
  { "displacement", topology::FieldBase::VECTOR, 3 },
  { "stress", topology::FieldBase::TENSOR, 6 },
  { "other", topology::FieldBase::OTHER, 2 },
};
const PylithScalar pylith::meshio::DataWriterHDF5DataFaultMeshHex8::_vertexFieldScalar[4*1] = {
  2.1, 3.2, 4.3, 5.4,
};
const PylithScalar pylith::meshio::DataWriterHDF5DataFaultMeshHex8::_vertexFieldVector[4*3] = {
  1.1, 2.2, 3.3,
  4.4, 5.5, 6.6,
  7.7, 8.8, 9.9,
  10.1, 11.2, 12.3,
};
const PylithScalar pylith::meshio::DataWriterHDF5DataFaultMeshHex8::_vertexFieldTensor[4*6] = {
  1.1, 1.2, 1.3, 1.4, 1.5, 1.6,
  2.1, 2.2, 2.3, 2.4, 2.5, 2.6,
  3.1, 3.2, 3.3, 3.4, 3.5, 3.6,
  4.1, 4.2, 4.3, 4.4, 4.5, 4.6,
};
const PylithScalar pylith::meshio::DataWriterHDF5DataFaultMeshHex8::_vertexFieldOther[4*2] = {
  1.2, 2.3,
  3.4, 4.5,
  5.6, 6.7,
  7.8, 8.9,
};

const int pylith::meshio::DataWriterHDF5DataFaultMeshHex8::_numCells = 1;

const pylith::meshio::DataWriterData::FieldStruct
pylith::meshio::DataWriterHDF5DataFaultMeshHex8::_cellFields[] = {
  { "pressure", topology::FieldBase::SCALAR, 1 },
  { "traction", topology::FieldBase::VECTOR, 3 },
  { "stress", topology::FieldBase::TENSOR, 6 },
  { "other", topology::FieldBase::OTHER, 3 },
};
const PylithScalar pylith::meshio::DataWriterHDF5DataFaultMeshHex8::_cellFieldScalar[1*1] = {
  2.1,
};
const PylithScalar pylith::meshio::DataWriterHDF5DataFaultMeshHex8::_cellFieldVector[1*3] = {
  1.1, 2.2, 3.3,
};
const PylithScalar pylith::meshio::DataWriterHDF5DataFaultMeshHex8::_cellFieldTensor[1*6] = {
  1.2, 2.3, 3.4, 4.5, 5.6, 6.7,
};
const PylithScalar pylith::meshio::DataWriterHDF5DataFaultMeshHex8::_cellFieldOther[1*3] = {
  1.2, 2.3, 3.4,
};

pylith::meshio::DataWriterHDF5DataFaultMeshHex8::DataWriterHDF5DataFaultMeshHex8(void)
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

pylith::meshio::DataWriterHDF5DataFaultMeshHex8::~DataWriterHDF5DataFaultMeshHex8(void)
{}


// End of file
