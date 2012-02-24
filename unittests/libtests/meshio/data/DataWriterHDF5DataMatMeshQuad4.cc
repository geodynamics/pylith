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

#include "DataWriterHDF5DataMatMeshQuad4.hh"

#include <assert.h> // USES assert()

const char* pylith::meshio::DataWriterHDF5DataMatMeshQuad4::_meshFilename = 
  "data/quad4.mesh";

const char* pylith::meshio::DataWriterHDF5DataMatMeshQuad4::_cellsLabel = 
  "material-id";
const int pylith::meshio::DataWriterHDF5DataMatMeshQuad4::_labelId = 2;

const char* pylith::meshio::DataWriterHDF5DataMatMeshQuad4::_timestepFilename = 
  "quad4_mat.h5";

const char* pylith::meshio::DataWriterHDF5DataMatMeshQuad4::_vertexFilename = 
  "quad4_mat_vertex.h5";

const char* pylith::meshio::DataWriterHDF5DataMatMeshQuad4::_cellFilename = 
  "quad4_mat_cell.h5";

const double pylith::meshio::DataWriterHDF5DataMatMeshQuad4::_time = 1.0;

const char* pylith::meshio::DataWriterHDF5DataMatMeshQuad4::_timeFormat = 
  "%3.1f";

const int pylith::meshio::DataWriterHDF5DataMatMeshQuad4::_numVertexFields = 3;
const int pylith::meshio::DataWriterHDF5DataMatMeshQuad4::_numVertices = 6;

const pylith::meshio::DataWriterData::FieldStruct
pylith::meshio::DataWriterHDF5DataMatMeshQuad4::_vertexFields[] = {
  { "displacements", topology::FieldBase::VECTOR, 2 },
  { "pressure", topology::FieldBase::SCALAR, 1 },
  { "other", topology::FieldBase::OTHER, 2 },
};
const double pylith::meshio::DataWriterHDF5DataMatMeshQuad4::_vertexField0[] = {
  1.1, 2.2,
  3.3, 4.4,
  5.5, 6.6,
  7.7, 8.8,
  9.9, 10.1,
  11.2, 12.3,
};
const double pylith::meshio::DataWriterHDF5DataMatMeshQuad4::_vertexField1[] = {
  2.1, 3.2, 4.3, 5.4, 6.5, 7.6
};
const double pylith::meshio::DataWriterHDF5DataMatMeshQuad4::_vertexField2[] = {
  1.2, 2.3,
  3.4, 4.5,
  5.6, 6.7,
  7.8, 8.9,
  9.8, 7.6,
  6.5, 5.4
};

const int pylith::meshio::DataWriterHDF5DataMatMeshQuad4::_numCellFields = 3;
const int pylith::meshio::DataWriterHDF5DataMatMeshQuad4::_numCells = 1;

const pylith::meshio::DataWriterData::FieldStruct
pylith::meshio::DataWriterHDF5DataMatMeshQuad4::_cellFields[] = {
  { "traction", topology::FieldBase::VECTOR, 2 },
  { "pressure", topology::FieldBase::SCALAR, 1 },
  { "other", topology::FieldBase::TENSOR, 3 },
};
const double pylith::meshio::DataWriterHDF5DataMatMeshQuad4::_cellField0[] = {
  1.1, 2.2,
};
const double pylith::meshio::DataWriterHDF5DataMatMeshQuad4::_cellField1[] = {
  2.1,
};
const double pylith::meshio::DataWriterHDF5DataMatMeshQuad4::_cellField2[] = {
  1.2, 2.3, 3.4,
};

pylith::meshio::DataWriterHDF5DataMatMeshQuad4::DataWriterHDF5DataMatMeshQuad4(void)
{ // constructor
  meshFilename = const_cast<char*>(_meshFilename);
  cellsLabel = const_cast<char*>(_cellsLabel);
  labelId = _labelId;

  timestepFilename = const_cast<char*>(_timestepFilename);
  vertexFilename = const_cast<char*>(_vertexFilename);
  cellFilename = const_cast<char*>(_cellFilename);

  time = _time;
  timeFormat = const_cast<char*>(_timeFormat);
  
  numVertexFields = _numVertexFields;
  assert(3 == numVertexFields);
  numVertices = _numVertices;
  vertexFieldsInfo = const_cast<DataWriterData::FieldStruct*>(_vertexFields);
  vertexFields[0] = const_cast<double*>(_vertexField0);
  vertexFields[1] = const_cast<double*>(_vertexField1);
  vertexFields[2] = const_cast<double*>(_vertexField2);

  numCellFields = _numCellFields;
  assert(3 == numCellFields);
  numCells = _numCells;
  cellFieldsInfo = const_cast<DataWriterData::FieldStruct*>(_cellFields);
  cellFields[0] = const_cast<double*>(_cellField0);
  cellFields[1] = const_cast<double*>(_cellField1);
  cellFields[2] = const_cast<double*>(_cellField2);
} // constructor

pylith::meshio::DataWriterHDF5DataMatMeshQuad4::~DataWriterHDF5DataMatMeshQuad4(void)
{}


// End of file
