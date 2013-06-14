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

#include "DataWriterHDF5DataMatMeshLine2.hh"

#include <assert.h> // USES assert()

const char* pylith::meshio::DataWriterHDF5DataMatMeshLine2::_meshFilename = 
  "data/line2.mesh";

const char* pylith::meshio::DataWriterHDF5DataMatMeshLine2::_cellsLabel = 
  "material-id";
const int pylith::meshio::DataWriterHDF5DataMatMeshLine2::_labelId = 0;

const char* pylith::meshio::DataWriterHDF5DataMatMeshLine2::_faultLabel = 
  "fault";
const int pylith::meshio::DataWriterHDF5DataMatMeshLine2::_faultId = 100;

const char* pylith::meshio::DataWriterHDF5DataMatMeshLine2::_timestepFilename = 
  "line2_mat.h5";

const char* pylith::meshio::DataWriterHDF5DataMatMeshLine2::_vertexFilename = 
  "line2_mat_vertex.h5";

const char* pylith::meshio::DataWriterHDF5DataMatMeshLine2::_cellFilename = 
  "line2_mat_cell.h5";

const PylithScalar pylith::meshio::DataWriterHDF5DataMatMeshLine2::_time = 1.0;

const char* pylith::meshio::DataWriterHDF5DataMatMeshLine2::_timeFormat = 
  "%3.1f";

const int pylith::meshio::DataWriterHDF5DataMatMeshLine2::_numVertices = 5;

const pylith::meshio::DataWriterData::FieldStruct
pylith::meshio::DataWriterHDF5DataMatMeshLine2::_vertexFields[4] = {
  { "pressure", topology::FieldBase::SCALAR, 1 },
  { "displacement", topology::FieldBase::VECTOR, 1 },
  { "stress", topology::FieldBase::TENSOR, 1 },
  { "other", topology::FieldBase::OTHER, 2 },
};
const PylithScalar pylith::meshio::DataWriterHDF5DataMatMeshLine2::_vertexFieldScalar[5*1] = {
  1.1, 2.2, 3.3, 4.4, 5.5
};
const PylithScalar pylith::meshio::DataWriterHDF5DataMatMeshLine2::_vertexFieldVector[5*1] = {
  2.1, 3.2, 4.3, 5.4, 6.5
};
const PylithScalar pylith::meshio::DataWriterHDF5DataMatMeshLine2::_vertexFieldTensor[5*1] = {
  3.1, 4.2, 5.3, 6.4, 7.5
};
const PylithScalar pylith::meshio::DataWriterHDF5DataMatMeshLine2::_vertexFieldOther[5*2] = {
  1.2, 2.3, 
  3.4, 4.5,
  5.6, 6.7,
  7.8, 8.9,
  9.1, 10.2
};

const int pylith::meshio::DataWriterHDF5DataMatMeshLine2::_numCells = 1;

const pylith::meshio::DataWriterData::FieldStruct
pylith::meshio::DataWriterHDF5DataMatMeshLine2::_cellFields[4] = {
  { "pressure", topology::FieldBase::SCALAR, 1 },
  { "traction", topology::FieldBase::VECTOR, 1 },
  { "stress", topology::FieldBase::TENSOR, 1 },
  { "other", topology::FieldBase::OTHER, 2 },
};
const PylithScalar pylith::meshio::DataWriterHDF5DataMatMeshLine2::_cellFieldScalar[1*1] = {
  1.1,
};
const PylithScalar pylith::meshio::DataWriterHDF5DataMatMeshLine2::_cellFieldVector[1*1] = {
  2.1,
};
const PylithScalar pylith::meshio::DataWriterHDF5DataMatMeshLine2::_cellFieldTensor[1*1] = {
  1.2,
};
const PylithScalar pylith::meshio::DataWriterHDF5DataMatMeshLine2::_cellFieldOther[1*2] = {
  1.2, 2.3,
};

pylith::meshio::DataWriterHDF5DataMatMeshLine2::DataWriterHDF5DataMatMeshLine2(void)
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

pylith::meshio::DataWriterHDF5DataMatMeshLine2::~DataWriterHDF5DataMatMeshLine2(void)
{}


// End of file
