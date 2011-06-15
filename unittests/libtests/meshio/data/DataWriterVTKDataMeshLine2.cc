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

#include "DataWriterVTKDataMeshLine2.hh"

#include <assert.h> // USES assert()

const char* pylith::meshio::DataWriterVTKDataMeshLine2::_meshFilename = 
  "data/line2.mesh";

const char* pylith::meshio::DataWriterVTKDataMeshLine2::_faultLabel = 
  "fault";
const int pylith::meshio::DataWriterVTKDataMeshLine2::_faultId = 100;

const char* pylith::meshio::DataWriterVTKDataMeshLine2::_timestepFilename = 
  "line2.vtk";

const char* pylith::meshio::DataWriterVTKDataMeshLine2::_vertexFilename = 
  "line2_vertex.vtk";

const char* pylith::meshio::DataWriterVTKDataMeshLine2::_cellFilename = 
  "line2_cell.vtk";

const double pylith::meshio::DataWriterVTKDataMeshLine2::_time = 1.0;

const char* pylith::meshio::DataWriterVTKDataMeshLine2::_timeFormat = 
  "%3.1f";

const int pylith::meshio::DataWriterVTKDataMeshLine2::_numVertexFields = 3;
const int pylith::meshio::DataWriterVTKDataMeshLine2::_numVertices = 5;

const pylith::meshio::DataWriterData::FieldStruct
pylith::meshio::DataWriterVTKDataMeshLine2::_vertexFields[] = {
  { "displacements", topology::FieldBase::VECTOR, 1 },
  { "pressure", topology::FieldBase::SCALAR, 1 },
  { "other", topology::FieldBase::OTHER, 2 },
};
const double pylith::meshio::DataWriterVTKDataMeshLine2::_vertexField0[] = {
  1.1, 2.2, 3.3, 4.4, 5.5
};
const double pylith::meshio::DataWriterVTKDataMeshLine2::_vertexField1[] = {
  2.1, 3.2, 4.3, 5.4, 6.5
};
const double pylith::meshio::DataWriterVTKDataMeshLine2::_vertexField2[] = {
  1.2, 2.3,
  3.4, 4.5,
  5.6, 6.7,
  7.8, 8.9,
  9.0, 10.1
};

const int pylith::meshio::DataWriterVTKDataMeshLine2::_numCellFields = 3;
const int pylith::meshio::DataWriterVTKDataMeshLine2::_numCells = 3;

const pylith::meshio::DataWriterData::FieldStruct
pylith::meshio::DataWriterVTKDataMeshLine2::_cellFields[] = {
  { "traction", topology::FieldBase::VECTOR, 1 },
  { "pressure", topology::FieldBase::SCALAR, 1 },
  { "other", topology::FieldBase::TENSOR, 1 },
};
const double pylith::meshio::DataWriterVTKDataMeshLine2::_cellField0[] = {
  1.1, 2.2, 3.3
};
const double pylith::meshio::DataWriterVTKDataMeshLine2::_cellField1[] = {
  2.1, 2.2, 2.3
};
const double pylith::meshio::DataWriterVTKDataMeshLine2::_cellField2[] = {
  1.2, 2.3, 3.4
};

pylith::meshio::DataWriterVTKDataMeshLine2::DataWriterVTKDataMeshLine2(void)
{ // constructor
  meshFilename = const_cast<char*>(_meshFilename);
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

pylith::meshio::DataWriterVTKDataMeshLine2::~DataWriterVTKDataMeshLine2(void)
{}


// End of file
