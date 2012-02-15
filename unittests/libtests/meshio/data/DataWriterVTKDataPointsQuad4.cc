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

#include "DataWriterVTKDataPointsQuad4.hh"

#include <assert.h> // USES assert()

const char* pylith::meshio::DataWriterVTKDataPointsQuad4::_meshFilename = 
  "data/quad4.mesh";

const char* pylith::meshio::DataWriterVTKDataPointsQuad4::_timestepFilename = 
  "quad4_points.vtk";

const char* pylith::meshio::DataWriterVTKDataPointsQuad4::_vertexFilename = 
  "quad4_points_vertex.vtk";

const PylithScalar pylith::meshio::DataWriterVTKDataPointsQuad4::_time = 1.0;

const char* pylith::meshio::DataWriterVTKDataPointsQuad4::_timeFormat = 
  "%3.1f";

const int pylith::meshio::DataWriterVTKDataPointsQuad4::_numVertexFields = 3;
const int pylith::meshio::DataWriterVTKDataPointsQuad4::_numVertices = 6;

const pylith::meshio::DataWriterData::FieldStruct
pylith::meshio::DataWriterVTKDataPointsQuad4::_vertexFields[] = {
  { "displacements", topology::FieldBase::VECTOR, 2 },
  { "pressure", topology::FieldBase::SCALAR, 1 },
  { "other", topology::FieldBase::OTHER, 2 },
};
const PylithScalar pylith::meshio::DataWriterVTKDataPointsQuad4::_vertexField0[] = {
  1.1, 2.2, // 0
  3.3, 4.4, // 1
  5.5, 6.6, // 2
  7.7, 8.8, // 3
  9.9, 10.1, // 4
  11.2, 12.3, // 5
};
const PylithScalar pylith::meshio::DataWriterVTKDataPointsQuad4::_vertexField1[] = {
  2.1, // 0
  3.2, // 1
  4.3, // 2
  5.4, // 3
  6.5, // 4
  7.6, // 5
};
const PylithScalar pylith::meshio::DataWriterVTKDataPointsQuad4::_vertexField2[] = {
  1.2, 2.3, // 0
  3.4, 4.5, // 1
  5.6, 6.7, // 2
  7.8, 8.9, // 3
  9.8, 7.6, // 4
  6.5, 5.4, // 5
};

const int pylith::meshio::DataWriterVTKDataPointsQuad4::_numPoints = 3;
const int pylith::meshio::DataWriterVTKDataPointsQuad4::_spaceDim = 2;
const PylithScalar pylith::meshio::DataWriterVTKDataPointsQuad4::_points[] = {
 -0.5, 0.0,
  0.00000001, 0.0,
  0.99999999, -0.99999999,
};


pylith::meshio::DataWriterVTKDataPointsQuad4::DataWriterVTKDataPointsQuad4(void)
{ // constructor
  meshFilename = const_cast<char*>(_meshFilename);

  timestepFilename = const_cast<char*>(_timestepFilename);
  vertexFilename = const_cast<char*>(_vertexFilename);

  time = _time;
  timeFormat = const_cast<char*>(_timeFormat);
  
  numVertexFields = _numVertexFields;
  numVertices = _numVertices;
  assert(3 == numVertexFields);
  vertexFieldsInfo = const_cast<DataWriterData::FieldStruct*>(_vertexFields);
  vertexFields[0] = const_cast<PylithScalar*>(_vertexField0);
  vertexFields[1] = const_cast<PylithScalar*>(_vertexField1);
  vertexFields[2] = const_cast<PylithScalar*>(_vertexField2);

  numPoints = _numPoints;
  spaceDim = _spaceDim;
  points = const_cast<PylithScalar*>(_points);
} // constructor


pylith::meshio::DataWriterVTKDataPointsQuad4::~DataWriterVTKDataPointsQuad4(void)
{}


// End of file
