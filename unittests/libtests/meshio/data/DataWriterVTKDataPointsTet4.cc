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

#include "DataWriterVTKDataPointsTet4.hh"

#include <assert.h> // USES assert()

const char* pylith::meshio::DataWriterVTKDataPointsTet4::_meshFilename = 
  "data/tet4.mesh";

const char* pylith::meshio::DataWriterVTKDataPointsTet4::_faultLabel = 
  "fault";
const int pylith::meshio::DataWriterVTKDataPointsTet4::_faultId = 100;

const char* pylith::meshio::DataWriterVTKDataPointsTet4::_timestepFilename = 
  "tet4_points.vtk";

const char* pylith::meshio::DataWriterVTKDataPointsTet4::_vertexFilename = 
  "tet4_points_vertex.vtk";

const PylithScalar pylith::meshio::DataWriterVTKDataPointsTet4::_time = 1.0;

const char* pylith::meshio::DataWriterVTKDataPointsTet4::_timeFormat = 
  "%3.1f";

const int pylith::meshio::DataWriterVTKDataPointsTet4::_numVertexFields = 3;
const int pylith::meshio::DataWriterVTKDataPointsTet4::_numVertices = 11;

const pylith::meshio::DataWriterData::FieldStruct
pylith::meshio::DataWriterVTKDataPointsTet4::_vertexFields[] = {
  { "displacements", topology::FieldBase::VECTOR, 3 },
  { "pressure", topology::FieldBase::SCALAR, 1 },
  { "other", topology::FieldBase::OTHER, 2 },
};
const PylithScalar pylith::meshio::DataWriterVTKDataPointsTet4::_vertexField0[] = {
  1.1, 2.2, 3.3, // 0
  4.4, 5.5, 6.6, // 1
  7.7, 8.8, 9.9, // 2
  10.0, 11.1, 12.2, // 3
  13.3, 14.4, 15.5, // 4
  16.6, 17.7, 18.8, // 5
  19.9, 20.0, 21.1, // 6
  22.2, 23.3, 24.4, // 7
  25.5, 26.6, 27.7, // 8
  28.8, 29.9, 30.0, // 9
  31.1, 32.2, 33.3, // 10
};
const PylithScalar pylith::meshio::DataWriterVTKDataPointsTet4::_vertexField1[] = {
  2.1, // 0
  3.2, // 1
  4.3, // 2
  5.4, // 3
  6.5, // 4
  7.6, // 5
  8.7, // 6
  9.8,  // 7
  10.9, // 8
  11.0, // 9
  12.1, // 10
};
const PylithScalar pylith::meshio::DataWriterVTKDataPointsTet4::_vertexField2[] = {
  1.2, 2.3, // 0
  3.4, 4.5, // 1
  5.6, 6.7, // 2
  7.8, 8.9, // 3
  9.0, 10.1, // 4
  11.2, 12.3, // 5
  13.4, 14.5, // 6
  15.6, 16.7, // 7
  17.8, 18.9, // 8
  19.0, 20.1, // 9
  21.2, 22.3, // 10
};

const int pylith::meshio::DataWriterVTKDataPointsTet4::_numPoints = 4;
const int pylith::meshio::DataWriterVTKDataPointsTet4::_spaceDim = 3;
const PylithScalar pylith::meshio::DataWriterVTKDataPointsTet4::_points[] = {
  -0.33333333, 0.0, 0.33333333,
  +0.00000001, 0.0, 0.33333333,
  +0.00000001, 0.0, 0.00000001,
  0.0, -0.99999999, 0.00000001,
};


pylith::meshio::DataWriterVTKDataPointsTet4::DataWriterVTKDataPointsTet4(void)
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


pylith::meshio::DataWriterVTKDataPointsTet4::~DataWriterVTKDataPointsTet4(void)
{}


// End of file
