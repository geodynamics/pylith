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

#include "DataWriterVTKDataPointsHex8.hh"

#include <assert.h> // USES assert()

const char* pylith::meshio::DataWriterVTKDataPointsHex8::_meshFilename = 
  "data/hex8.mesh";

const char* pylith::meshio::DataWriterVTKDataPointsHex8::_faultLabel = 
  "fault";
const int pylith::meshio::DataWriterVTKDataPointsHex8::_faultId = 100;

const char* pylith::meshio::DataWriterVTKDataPointsHex8::_timestepFilename = 
  "hex8_points.vtk";

const char* pylith::meshio::DataWriterVTKDataPointsHex8::_vertexFilename = 
  "hex8_points_vertex.vtk";

const PylithScalar pylith::meshio::DataWriterVTKDataPointsHex8::_time = 1.0;

const char* pylith::meshio::DataWriterVTKDataPointsHex8::_timeFormat = 
  "%3.1f";

const int pylith::meshio::DataWriterVTKDataPointsHex8::_numVertexFields = 3;
const int pylith::meshio::DataWriterVTKDataPointsHex8::_numVertices = 20;

const pylith::meshio::DataWriterData::FieldStruct
pylith::meshio::DataWriterVTKDataPointsHex8::_vertexFields[] = {
  { "displacements", topology::FieldBase::VECTOR, 2 },
  { "pressure", topology::FieldBase::SCALAR, 1 },
  { "other", topology::FieldBase::OTHER, 2 },
};
const PylithScalar pylith::meshio::DataWriterVTKDataPointsHex8::_vertexField0[] = {
  1.1, 2.2, 3.3, // 0
  4.4, 5.5, 6.6, // 1
  7.7, 8.8, 9.9, // 2
  10.1, 11.2, 12.3, // 3
  1.2, 2.3, 3.4, // 4
  4.5, 5.6, 6.7, // 5
  7.8, 8.9, 9.0, // 6
  10.2, 11.3, 12.4, // 7
  1.3, 2.4, 3.5, // 8
  4.6, 5.7, 6.8, // 9
  7.9, 8.0, 9.1, // 10
  10.2, 11.3, 12.4, // 11
  13.5, 14.6, 15.7, // 12
  16.8, 17.9, 18.1, // 13
  19.2, 20.3, 21.4, // 14
  22.5, 23.6, 24.7, // 15
  25.8, 26.9, 27.1, // 16
  28.8, 29.9, 30.1, // 17
  31.8, 32.9, 33.1, // 18
  34.8, 35.9, 36.1, // 19
};
const PylithScalar pylith::meshio::DataWriterVTKDataPointsHex8::_vertexField1[] = {
  2.1, // 0
  3.2, // 1
  4.3, // 2 
  5.4, // 3
  6.5, // 4
  7.6, // 5
  8.7, // 6
  9.8, // 7
  10.0, // 8
  12.1, // 9
  11.1, // 10
  13.1, // 11
  14.1, // 12
  15.1, // 13
  16.1, // 14
  17.1, // 15
  18.1, // 16
  19.1, // 17
  20.1, // 18
  21.2, // 19
};
const PylithScalar pylith::meshio::DataWriterVTKDataPointsHex8::_vertexField2[] = {
  1.2, 2.3, // 0
  3.4, 4.5, // 1
  5.6, 6.7, // 2
  7.8, 8.9, // 3
  1.3, 2.4, // 4
  3.5, 4.6, // 5
  5.7, 6.8, // 6
  7.9, 8.0, // 7
  1.3, 2.4, // 8
  3.5, 4.6, // 9
  5.7, 6.8, // 10
  8.0, 1.4, // 11
  2.5, 3.6, // 12
  4.8, 1.5, // 13
  2.6, 3.7, // 14
  4.8, 5.9, // 15
  6.1, 7.2, // 16
  7.1, 8.2, // 17
  8.1, 9.2, // 18
  9.1, 10.1, // 19
};

const int pylith::meshio::DataWriterVTKDataPointsHex8::_numPoints = 4;
const int pylith::meshio::DataWriterVTKDataPointsHex8::_spaceDim = 3;
const PylithScalar pylith::meshio::DataWriterVTKDataPointsHex8::_points[] = {
  -0.5, 0.0, 0.5,
  -0.00000001, 0.0, 0.0,
  -0.00000001, 0.0, 0.99999999,
   0.99999999, 0.99999999, -0.99999999,
};


pylith::meshio::DataWriterVTKDataPointsHex8::DataWriterVTKDataPointsHex8(void)
{ // constructor
  meshFilename = const_cast<char*>(_meshFilename);
  faultLabel = const_cast<char*>(_faultLabel);
  faultId = _faultId;

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


pylith::meshio::DataWriterVTKDataPointsHex8::~DataWriterVTKDataPointsHex8(void)
{}


// End of file
