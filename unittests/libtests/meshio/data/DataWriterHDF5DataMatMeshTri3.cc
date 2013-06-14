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

#include "DataWriterHDF5DataMatMeshTri3.hh"

#include <assert.h> // USES assert()

const char* pylith::meshio::DataWriterHDF5DataMatMeshTri3::_meshFilename = 
  "data/tri3.mesh";

const char* pylith::meshio::DataWriterHDF5DataMatMeshTri3::_cellsLabel = 
  "material-id";
const int pylith::meshio::DataWriterHDF5DataMatMeshTri3::_labelId = 0;

const char* pylith::meshio::DataWriterHDF5DataMatMeshTri3::_faultLabel = 
  "fault";
const int pylith::meshio::DataWriterHDF5DataMatMeshTri3::_faultId = 100;

const char* pylith::meshio::DataWriterHDF5DataMatMeshTri3::_timestepFilename = 
  "tri3_mat.h5";

const char* pylith::meshio::DataWriterHDF5DataMatMeshTri3::_vertexFilename = 
  "tri3_mat_vertex.h5";

const char* pylith::meshio::DataWriterHDF5DataMatMeshTri3::_cellFilename = 
  "tri3_mat_cell.h5";

const PylithScalar pylith::meshio::DataWriterHDF5DataMatMeshTri3::_time = 1.0;

const char* pylith::meshio::DataWriterHDF5DataMatMeshTri3::_timeFormat = 
  "%3.1f";

const int pylith::meshio::DataWriterHDF5DataMatMeshTri3::_numVertices = 8;

const pylith::meshio::DataWriterData::FieldStruct
pylith::meshio::DataWriterHDF5DataMatMeshTri3::_vertexFields[4] = {
  { "pressure", topology::FieldBase::SCALAR, 1 },
  { "displacement", topology::FieldBase::VECTOR, 2 },
  { "stress", topology::FieldBase::TENSOR, 3 },
  { "other", topology::FieldBase::OTHER, 2 },
};
const PylithScalar pylith::meshio::DataWriterHDF5DataMatMeshTri3::_vertexFieldVector[8*2] = {
  1.1, 2.2,
  3.3, 4.4,
  5.5, 6.6,
  7.7, 8.8,
  9.9, 10.0,
  11.1, 12.2,
  13.3, 14.4,
  15.5, 16.6,
};
const PylithScalar pylith::meshio::DataWriterHDF5DataMatMeshTri3::_vertexFieldScalar[8*1] = {
  2.1, 3.2, 4.3, 5.4, 6.5, 7.6, 8.7, 9.8
};
const PylithScalar pylith::meshio::DataWriterHDF5DataMatMeshTri3::_vertexFieldTensor[8*3] = {
  1.1, 1.2, 1.3,
  2.1, 2.2, 3.3,
  3.1, 3.2, 4.3,
  4.1, 4.2, 5.3,
  5.1, 5.2, 6.3,
  6.1, 6.2, 7.3,
  7.1, 7.2, 8.3,
  8.1, 8.2, 9.3,
};
const PylithScalar pylith::meshio::DataWriterHDF5DataMatMeshTri3::_vertexFieldOther[8*2] = {
  1.2, 2.3,
  3.4, 4.5,
  5.6, 6.7,
  7.8, 8.9,
  9.0, 10.1,
  11.2, 12.3,
  13.4, 14.5,
  15.6, 16.7
};

const int pylith::meshio::DataWriterHDF5DataMatMeshTri3::_numCells = 1;

const pylith::meshio::DataWriterData::FieldStruct
pylith::meshio::DataWriterHDF5DataMatMeshTri3::_cellFields[4] = {
  { "pressure", topology::FieldBase::SCALAR, 1 },
  { "traction", topology::FieldBase::VECTOR, 2 },
  { "stress", topology::FieldBase::TENSOR, 3 },
  { "other", topology::FieldBase::OTHER, 2 },
};
const PylithScalar pylith::meshio::DataWriterHDF5DataMatMeshTri3::_cellFieldScalar[1*1] = {
  2.1,
};
const PylithScalar pylith::meshio::DataWriterHDF5DataMatMeshTri3::_cellFieldVector[1*2] = {
  1.1, 2.2,
};
const PylithScalar pylith::meshio::DataWriterHDF5DataMatMeshTri3::_cellFieldTensor[1*3] = {
  1.2, 2.3, 3.4,
};
const PylithScalar pylith::meshio::DataWriterHDF5DataMatMeshTri3::_cellFieldOther[1*2] = {
  1.2, 2.3,
};

pylith::meshio::DataWriterHDF5DataMatMeshTri3::DataWriterHDF5DataMatMeshTri3(void)
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

pylith::meshio::DataWriterHDF5DataMatMeshTri3::~DataWriterHDF5DataMatMeshTri3(void)
{}


// End of file
