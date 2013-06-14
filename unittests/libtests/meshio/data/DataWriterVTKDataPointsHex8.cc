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

const int pylith::meshio::DataWriterVTKDataPointsHex8::_numVertices = 20;

const pylith::meshio::DataWriterData::FieldStruct
pylith::meshio::DataWriterVTKDataPointsHex8::_vertexFields[4] = {
  { "pressure", topology::FieldBase::SCALAR, 1 },
  { "displacement", topology::FieldBase::VECTOR, 3 },
  { "stress", topology::FieldBase::TENSOR, 6 },
  { "other", topology::FieldBase::OTHER, 2 },
};
const PylithScalar pylith::meshio::DataWriterVTKDataPointsHex8::_vertexFieldScalar[20*1] = {
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
const PylithScalar pylith::meshio::DataWriterVTKDataPointsHex8::_vertexFieldVector[20*3] = {
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
const PylithScalar pylith::meshio::DataWriterVTKDataPointsHex8::_vertexFieldTensor[20*6] = {
  1.1, 1.2, 1.3, 1.4, 1.5, 1.6,
  2.1, 2.2, 2.3, 2.4, 2.5, 2.6,
  3.1, 3.2, 3.3, 3.4, 3.5, 3.6,
  4.1, 4.2, 4.3, 4.4, 4.5, 4.6,
  5.1, 5.2, 5.3, 5.4, 5.5, 5.6,
  6.1, 6.2, 6.3, 6.4, 6.5, 6.6,
  7.1, 7.2, 7.3, 7.4, 7.5, 7.6,
  8.1, 8.2, 8.3, 8.4, 8.5, 8.6,
  9.1, 9.2, 9.3, 9.4, 9.5, 9.6,
  10.1, 10.2, 10.3, 10.4, 10.5, 10.6,
  11.1, 11.2, 11.3, 11.4, 11.5, 11.6,
  12.1, 12.2, 12.3, 12.4, 12.5, 12.6,
  13.1, 13.2, 13.3, 13.4, 13.5, 13.6,
  14.1, 14.2, 14.3, 14.4, 14.5, 14.6,
  15.1, 15.2, 15.3, 15.4, 15.5, 15.6,
  16.1, 16.2, 16.3, 16.4, 16.5, 16.6,
  17.1, 17.2, 17.3, 17.4, 17.5, 17.6,
  18.1, 18.2, 18.3, 18.4, 18.5, 18.6,
  19.1, 19.2, 19.3, 19.4, 19.5, 19.6,
  20.1, 20.2, 20.3, 20.4, 20.5, 20.6,
};
const PylithScalar pylith::meshio::DataWriterVTKDataPointsHex8::_vertexFieldOther[20*2] = {
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
const PylithScalar pylith::meshio::DataWriterVTKDataPointsHex8::_points[4*3] = {
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
  
  numVertices = _numVertices;
  assert(DataWriterData::numVertexFields == numVertexFields);
  vertexFieldsInfo = const_cast<DataWriterData::FieldStruct*>(_vertexFields);
  vertexFields[0] = const_cast<PylithScalar*>(_vertexFieldScalar);
  vertexFields[1] = const_cast<PylithScalar*>(_vertexFieldVector);
  vertexFields[2] = const_cast<PylithScalar*>(_vertexFieldTensor);
  vertexFields[3] = const_cast<PylithScalar*>(_vertexFieldOther);

  numPoints = _numPoints;
  spaceDim = _spaceDim;
  points = const_cast<PylithScalar*>(_points);
} // constructor


pylith::meshio::DataWriterVTKDataPointsHex8::~DataWriterVTKDataPointsHex8(void)
{}


// End of file
