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

#include "DataWriterVTKDataFaultMeshTri3.hh"

#include <assert.h> // USES assert()

const char* pylith::meshio::DataWriterVTKDataFaultMeshTri3::_meshFilename = 
  "data/tri3.mesh";

const char* pylith::meshio::DataWriterVTKDataFaultMeshTri3::_faultLabel = 
  "fault";
const int pylith::meshio::DataWriterVTKDataFaultMeshTri3::_faultId = 100;

const char* pylith::meshio::DataWriterVTKDataFaultMeshTri3::_timestepFilename = 
  "tri3_fault.vtk";

const char* pylith::meshio::DataWriterVTKDataFaultMeshTri3::_vertexFilename = 
  "tri3_fault_vertex.vtk";

const char* pylith::meshio::DataWriterVTKDataFaultMeshTri3::_cellFilename = 
  "tri3_fault_cell.vtk";

const PylithScalar pylith::meshio::DataWriterVTKDataFaultMeshTri3::_time = 1.0;

const char* pylith::meshio::DataWriterVTKDataFaultMeshTri3::_timeFormat = 
  "%3.1f";

const int pylith::meshio::DataWriterVTKDataFaultMeshTri3::_numVertexFields = 3;
const int pylith::meshio::DataWriterVTKDataFaultMeshTri3::_numVertices = 2;

const pylith::meshio::DataWriterData::FieldStruct
pylith::meshio::DataWriterVTKDataFaultMeshTri3::_vertexFields[] = {
  { "displacements", topology::FieldBase::VECTOR, 2 },
  { "pressure", topology::FieldBase::SCALAR, 1 },
  { "other", topology::FieldBase::OTHER, 2 },
};
const PylithScalar pylith::meshio::DataWriterVTKDataFaultMeshTri3::_vertexField0[] = {
  1.1, 2.2,
  3.3, 4.4,
};
const PylithScalar pylith::meshio::DataWriterVTKDataFaultMeshTri3::_vertexField1[] = {
  2.1, 3.2,
};
const PylithScalar pylith::meshio::DataWriterVTKDataFaultMeshTri3::_vertexField2[] = {
  1.2, 2.3,
  3.4, 4.5,
};

const int pylith::meshio::DataWriterVTKDataFaultMeshTri3::_numCellFields = 3;
const int pylith::meshio::DataWriterVTKDataFaultMeshTri3::_numCells = 1;

const pylith::meshio::DataWriterData::FieldStruct
pylith::meshio::DataWriterVTKDataFaultMeshTri3::_cellFields[] = {
  { "traction", topology::FieldBase::VECTOR, 2 },
  { "pressure", topology::FieldBase::SCALAR, 1 },
  { "other", topology::FieldBase::TENSOR, 3 },
};
const PylithScalar pylith::meshio::DataWriterVTKDataFaultMeshTri3::_cellField0[] = {
  1.1, 2.2,
};
const PylithScalar pylith::meshio::DataWriterVTKDataFaultMeshTri3::_cellField1[] = {
  2.1,
};
const PylithScalar pylith::meshio::DataWriterVTKDataFaultMeshTri3::_cellField2[] = {
  1.2, 2.3, 3.4,
};

pylith::meshio::DataWriterVTKDataFaultMeshTri3::DataWriterVTKDataFaultMeshTri3(void)
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

pylith::meshio::DataWriterVTKDataFaultMeshTri3::~DataWriterVTKDataFaultMeshTri3(void)
{}


// End of file
