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

#include "DataWriterHDF5DataMatMeshHex8.hh"

#include <assert.h> // USES assert()

const char* pylith::meshio::DataWriterHDF5DataMatMeshHex8::_meshFilename = 
  "data/hex8.mesh";

const char* pylith::meshio::DataWriterHDF5DataMatMeshHex8::_cellsLabel = 
  "material-id";
const int pylith::meshio::DataWriterHDF5DataMatMeshHex8::_labelId = 0;

const char* pylith::meshio::DataWriterHDF5DataMatMeshHex8::_faultLabel = 
  "fault";
const int pylith::meshio::DataWriterHDF5DataMatMeshHex8::_faultId = 100;

const char* pylith::meshio::DataWriterHDF5DataMatMeshHex8::_timestepFilename = 
  "hex8_mat.h5";

const char* pylith::meshio::DataWriterHDF5DataMatMeshHex8::_vertexFilename = 
  "hex8_mat_vertex.h5";

const char* pylith::meshio::DataWriterHDF5DataMatMeshHex8::_cellFilename = 
  "hex8_mat_cell.h5";

const double pylith::meshio::DataWriterHDF5DataMatMeshHex8::_time = 1.0;

const char* pylith::meshio::DataWriterHDF5DataMatMeshHex8::_timeFormat = 
  "%3.1f";

const int pylith::meshio::DataWriterHDF5DataMatMeshHex8::_numVertexFields = 3;
const int pylith::meshio::DataWriterHDF5DataMatMeshHex8::_numVertices = 20;

const pylith::meshio::DataWriterData::FieldStruct
pylith::meshio::DataWriterHDF5DataMatMeshHex8::_vertexFields[] = {
  { "displacements", topology::FieldBase::VECTOR, 3 },
  { "pressure", topology::FieldBase::SCALAR, 1 },
  { "other", topology::FieldBase::OTHER, 2 },
};
const double pylith::meshio::DataWriterHDF5DataMatMeshHex8::_vertexField0[] = {
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
  25.5, 26.6, 27.7,
  28.8, 29.9, 30.0,
  31.1, 32.2, 33.3,
  34.4, 35.5, 36.6
};
const double pylith::meshio::DataWriterHDF5DataMatMeshHex8::_vertexField1[] = {
  2.1, 3.2, 4.3, 5.4, 6.5, 7.6, 8.7, 9.8,
  10.0, 12.1, 11.1, 13.1, 14.1, 15.1, 16.1, 17.1,
  18.1, 19.1, 20.1, 21.1
};
const double pylith::meshio::DataWriterHDF5DataMatMeshHex8::_vertexField2[] = {
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
  1.1, 1.2,
  2.1, 2.2,
  3.1, 3.2,
  4.1, 4.2,
};

const int pylith::meshio::DataWriterHDF5DataMatMeshHex8::_numCellFields = 3;
const int pylith::meshio::DataWriterHDF5DataMatMeshHex8::_numCells = 1;

const pylith::meshio::DataWriterData::FieldStruct
pylith::meshio::DataWriterHDF5DataMatMeshHex8::_cellFields[] = {
  { "traction", topology::FieldBase::VECTOR, 3 },
  { "pressure", topology::FieldBase::SCALAR, 1 },
  { "other", topology::FieldBase::TENSOR, 6 },
};
const double pylith::meshio::DataWriterHDF5DataMatMeshHex8::_cellField0[] = {
  1.1, 2.2, 3.3,
};
const double pylith::meshio::DataWriterHDF5DataMatMeshHex8::_cellField1[] = {
  2.1,
};
const double pylith::meshio::DataWriterHDF5DataMatMeshHex8::_cellField2[] = {
  1.2, 2.3, 3.4, 4.5, 5.6, 6.7,
};

pylith::meshio::DataWriterHDF5DataMatMeshHex8::DataWriterHDF5DataMatMeshHex8(void)
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
  
  numVertexFields = _numVertexFields;
  numVertices = _numVertices;
  assert(3 == numVertexFields);
  vertexFieldsInfo = const_cast<DataWriterData::FieldStruct*>(_vertexFields);
  vertexFields[0] = const_cast<double*>(_vertexField0);
  vertexFields[1] = const_cast<double*>(_vertexField1);
  vertexFields[2] = const_cast<double*>(_vertexField2);

  numCellFields = _numCellFields;
  numCells = _numCells;
  assert(3 == numCellFields);
  cellFieldsInfo = const_cast<DataWriterData::FieldStruct*>(_cellFields);
  cellFields[0] = const_cast<double*>(_cellField0);
  cellFields[1] = const_cast<double*>(_cellField1);
  cellFields[2] = const_cast<double*>(_cellField2);
} // constructor

pylith::meshio::DataWriterHDF5DataMatMeshHex8::~DataWriterHDF5DataMatMeshHex8(void)
{}


// End of file
