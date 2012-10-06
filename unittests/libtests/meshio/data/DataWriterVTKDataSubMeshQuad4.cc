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

#include "DataWriterVTKDataSubMeshQuad4.hh"

#include <assert.h> // USES assert()

const char* pylith::meshio::DataWriterVTKDataSubMeshQuad4::_meshFilename = 
  "data/quad4.mesh";

const char* pylith::meshio::DataWriterVTKDataSubMeshQuad4::_bcLabel = 
  "bc3";

const char* pylith::meshio::DataWriterVTKDataSubMeshQuad4::_timestepFilename = 
  "quad4_surf.vtk";

const char* pylith::meshio::DataWriterVTKDataSubMeshQuad4::_vertexFilename = 
  "quad4_surf_vertex.vtk";

const char* pylith::meshio::DataWriterVTKDataSubMeshQuad4::_cellFilename = 
  "quad4_surf_cell.vtk";

const PylithScalar pylith::meshio::DataWriterVTKDataSubMeshQuad4::_time = 1.0;

const char* pylith::meshio::DataWriterVTKDataSubMeshQuad4::_timeFormat = 
  "%3.1f";

const int pylith::meshio::DataWriterVTKDataSubMeshQuad4::_numVertexFields = 3;
const int pylith::meshio::DataWriterVTKDataSubMeshQuad4::_numVertices = 3;

const pylith::meshio::DataWriterData::FieldStruct
pylith::meshio::DataWriterVTKDataSubMeshQuad4::_vertexFields[] = {
  { "displacements", topology::FieldBase::VECTOR, 2 },
  { "pressure", topology::FieldBase::SCALAR, 1 },
  { "other", topology::FieldBase::OTHER, 2 },
};
const PylithScalar pylith::meshio::DataWriterVTKDataSubMeshQuad4::_vertexField0[] = {
  1.1, 2.2,
  5.5, 6.6,
  9.9, 10.0,
};
const PylithScalar pylith::meshio::DataWriterVTKDataSubMeshQuad4::_vertexField1[] = {
  2.1, 4.3, 6.5,
};
const PylithScalar pylith::meshio::DataWriterVTKDataSubMeshQuad4::_vertexField2[] = {
  1.2, 2.3,
  5.6, 6.7,
  9.1, 10.2,
};

const int pylith::meshio::DataWriterVTKDataSubMeshQuad4::_numCellFields = 3;
const int pylith::meshio::DataWriterVTKDataSubMeshQuad4::_numCells = 2;

const pylith::meshio::DataWriterData::FieldStruct
pylith::meshio::DataWriterVTKDataSubMeshQuad4::_cellFields[] = {
  { "traction", topology::FieldBase::VECTOR, 2 },
  { "pressure", topology::FieldBase::SCALAR, 1 },
  { "other", topology::FieldBase::TENSOR, 3 },
};
const PylithScalar pylith::meshio::DataWriterVTKDataSubMeshQuad4::_cellField0[] = {
  1.1, 2.2,
  3.3, 4.4,
};
const PylithScalar pylith::meshio::DataWriterVTKDataSubMeshQuad4::_cellField1[] = {
  2.1, 3.2,
};
const PylithScalar pylith::meshio::DataWriterVTKDataSubMeshQuad4::_cellField2[] = {
  1.2, 2.3, 3.4,
  4.5, 5.6, 6.7,
};

pylith::meshio::DataWriterVTKDataSubMeshQuad4::DataWriterVTKDataSubMeshQuad4(void)
{ // constructor
  meshFilename = const_cast<char*>(_meshFilename);
  bcLabel = const_cast<char*>(_bcLabel);

  timestepFilename = const_cast<char*>(_timestepFilename);
  vertexFilename = const_cast<char*>(_vertexFilename);
  cellFilename = const_cast<char*>(_cellFilename);

  time = _time;
  timeFormat = const_cast<char*>(_timeFormat);
  
  numVertexFields = _numVertexFields;
  assert(3 == numVertexFields);
  numVertices = _numVertices;
  vertexFieldsInfo = const_cast<DataWriterData::FieldStruct*>(_vertexFields);
  vertexFields[0] = const_cast<PylithScalar*>(_vertexField0);
  vertexFields[1] = const_cast<PylithScalar*>(_vertexField1);
  vertexFields[2] = const_cast<PylithScalar*>(_vertexField2);

  numCellFields = _numCellFields;
  assert(3 == numCellFields);
  numCells = _numCells;
  cellFieldsInfo = const_cast<DataWriterData::FieldStruct*>(_cellFields);
  cellFields[0] = const_cast<PylithScalar*>(_cellField0);
  cellFields[1] = const_cast<PylithScalar*>(_cellField1);
  cellFields[2] = const_cast<PylithScalar*>(_cellField2);
} // constructor

pylith::meshio::DataWriterVTKDataSubMeshQuad4::~DataWriterVTKDataSubMeshQuad4(void)
{}


// End of file
