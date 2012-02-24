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

#include "DataWriterVTKDataMatMeshLine2.hh"

#include <assert.h> // USES assert()

const char* pylith::meshio::DataWriterVTKDataMatMeshLine2::_meshFilename = 
  "data/line2.mesh";

const char* pylith::meshio::DataWriterVTKDataMatMeshLine2::_cellsLabel = 
  "material-id";
const int pylith::meshio::DataWriterVTKDataMatMeshLine2::_labelId = 0;

const char* pylith::meshio::DataWriterVTKDataMatMeshLine2::_faultLabel = 
  "fault";
const int pylith::meshio::DataWriterVTKDataMatMeshLine2::_faultId = 100;

const char* pylith::meshio::DataWriterVTKDataMatMeshLine2::_timestepFilename = 
  "line2_mat.vtk";

const char* pylith::meshio::DataWriterVTKDataMatMeshLine2::_vertexFilename = 
  "line2_mat_vertex.vtk";

const char* pylith::meshio::DataWriterVTKDataMatMeshLine2::_cellFilename = 
  "line2_mat_cell.vtk";

const PylithScalar pylith::meshio::DataWriterVTKDataMatMeshLine2::_time = 1.0;

const char* pylith::meshio::DataWriterVTKDataMatMeshLine2::_timeFormat = 
  "%3.1f";

const int pylith::meshio::DataWriterVTKDataMatMeshLine2::_numVertexFields = 3;
const int pylith::meshio::DataWriterVTKDataMatMeshLine2::_numVertices = 5;

const pylith::meshio::DataWriterData::FieldStruct
pylith::meshio::DataWriterVTKDataMatMeshLine2::_vertexFields[] = {
  { "displacements", topology::FieldBase::VECTOR, 1 },
  { "pressure", topology::FieldBase::SCALAR, 1 },
  { "other", topology::FieldBase::OTHER, 2 },
};
const PylithScalar pylith::meshio::DataWriterVTKDataMatMeshLine2::_vertexField0[] = {
  1.1, 2.2, 3.3, 4.4, 5.5
};
const PylithScalar pylith::meshio::DataWriterVTKDataMatMeshLine2::_vertexField1[] = {
  2.1, 3.2, 4.3, 5.4, 6.5
};
const PylithScalar pylith::meshio::DataWriterVTKDataMatMeshLine2::_vertexField2[] = {
  1.2, 2.3, 
  3.4, 4.5,
  5.6, 6.7,
  7.8, 8.9,
  9.1, 10.2
};

const int pylith::meshio::DataWriterVTKDataMatMeshLine2::_numCellFields = 3;
const int pylith::meshio::DataWriterVTKDataMatMeshLine2::_numCells = 1;

const pylith::meshio::DataWriterData::FieldStruct
pylith::meshio::DataWriterVTKDataMatMeshLine2::_cellFields[] = {
  { "traction", topology::FieldBase::VECTOR, 1 },
  { "pressure", topology::FieldBase::SCALAR, 1 },
  { "other", topology::FieldBase::TENSOR, 1 },
};
const PylithScalar pylith::meshio::DataWriterVTKDataMatMeshLine2::_cellField0[] = {
  1.1,
};
const PylithScalar pylith::meshio::DataWriterVTKDataMatMeshLine2::_cellField1[] = {
  2.1,
};
const PylithScalar pylith::meshio::DataWriterVTKDataMatMeshLine2::_cellField2[] = {
  1.2,
};

pylith::meshio::DataWriterVTKDataMatMeshLine2::DataWriterVTKDataMatMeshLine2(void)
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
  vertexFields[0] = const_cast<PylithScalar*>(_vertexField0);
  vertexFields[1] = const_cast<PylithScalar*>(_vertexField1);
  vertexFields[2] = const_cast<PylithScalar*>(_vertexField2);

  numCellFields = _numCellFields;
  numCells = _numCells;
  assert(3 == numCellFields);
  cellFieldsInfo = const_cast<DataWriterData::FieldStruct*>(_cellFields);
  cellFields[0] = const_cast<PylithScalar*>(_cellField0);
  cellFields[1] = const_cast<PylithScalar*>(_cellField1);
  cellFields[2] = const_cast<PylithScalar*>(_cellField2);
} // constructor

pylith::meshio::DataWriterVTKDataMatMeshLine2::~DataWriterVTKDataMatMeshLine2(void)
{}


// End of file
