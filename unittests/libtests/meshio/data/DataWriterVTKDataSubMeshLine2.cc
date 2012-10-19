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

#include "DataWriterVTKDataSubMeshLine2.hh"

#include <assert.h> // USES assert()

const char* pylith::meshio::DataWriterVTKDataSubMeshLine2::_meshFilename = 
  "data/line2.mesh";

const char* pylith::meshio::DataWriterVTKDataSubMeshLine2::_bcLabel = 
  "bc1";

const char* pylith::meshio::DataWriterVTKDataSubMeshLine2::_faultLabel = 
  "fault";
const int pylith::meshio::DataWriterVTKDataSubMeshLine2::_faultId = 100;

const char* pylith::meshio::DataWriterVTKDataSubMeshLine2::_timestepFilename = 
  "line2_surf.vtk";

const char* pylith::meshio::DataWriterVTKDataSubMeshLine2::_vertexFilename = 
  "line2_surf_vertex.vtk";

const char* pylith::meshio::DataWriterVTKDataSubMeshLine2::_cellFilename = 
  "line2_surf_cell.vtk";

const PylithScalar pylith::meshio::DataWriterVTKDataSubMeshLine2::_time = 1.0;

const char* pylith::meshio::DataWriterVTKDataSubMeshLine2::_timeFormat = 
  "%3.1f";

const int pylith::meshio::DataWriterVTKDataSubMeshLine2::_numVertices = 1;

const pylith::meshio::DataWriterData::FieldStruct
pylith::meshio::DataWriterVTKDataSubMeshLine2::_vertexFields[4] = {
  { "pressure", topology::FieldBase::SCALAR, 1 },
  { "displacement", topology::FieldBase::VECTOR, 1 },
  { "stress", topology::FieldBase::TENSOR, 1 },
  { "other", topology::FieldBase::OTHER, 2 },
};
const PylithScalar pylith::meshio::DataWriterVTKDataSubMeshLine2::_vertexFieldScalar[5*1] = {
  1.1,  3.3,  5.5,  7.7,  9.9,
};
const PylithScalar pylith::meshio::DataWriterVTKDataSubMeshLine2::_vertexFieldVector[5*1] = {
  2.1
};
const PylithScalar pylith::meshio::DataWriterVTKDataSubMeshLine2::_vertexFieldTensor[5*1] = {
  3.1, 4.2, 5.3, 6.4, 7.5,
};
const PylithScalar pylith::meshio::DataWriterVTKDataSubMeshLine2::_vertexFieldOther[5*2] = {
  1.2, 2.3,
};

const int pylith::meshio::DataWriterVTKDataSubMeshLine2::_numCells = 1;

const pylith::meshio::DataWriterData::FieldStruct
pylith::meshio::DataWriterVTKDataSubMeshLine2::_cellFields[4] = {
  { "pressure", topology::FieldBase::SCALAR, 1 },
  { "traction", topology::FieldBase::VECTOR, 1 },
  { "stress", topology::FieldBase::TENSOR, 1 },
  { "other", topology::FieldBase::OTHER, 3 },
};
const PylithScalar pylith::meshio::DataWriterVTKDataSubMeshLine2::_cellFieldScalar[1*1] = {
  1.1,
};
const PylithScalar pylith::meshio::DataWriterVTKDataSubMeshLine2::_cellFieldVector[1*1] = {
  2.1,
};
const PylithScalar pylith::meshio::DataWriterVTKDataSubMeshLine2::_cellFieldTensor[1*1] = {
  3.1,
};
const PylithScalar pylith::meshio::DataWriterVTKDataSubMeshLine2::_cellFieldOther[1*3] = {
  1.2, 2.3, 3.4,
};

pylith::meshio::DataWriterVTKDataSubMeshLine2::DataWriterVTKDataSubMeshLine2(void)
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

pylith::meshio::DataWriterVTKDataSubMeshLine2::~DataWriterVTKDataSubMeshLine2(void)
{}


// End of file
