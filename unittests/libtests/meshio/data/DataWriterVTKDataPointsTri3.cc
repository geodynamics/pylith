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

#include "DataWriterVTKDataPointsTri3.hh"

#include <assert.h> // USES assert()

const char* pylith::meshio::DataWriterVTKDataPointsTri3::_meshFilename = 
  "data/tri3.mesh";

const char* pylith::meshio::DataWriterVTKDataPointsTri3::_faultLabel = 
  "fault";
const int pylith::meshio::DataWriterVTKDataPointsTri3::_faultId = 100;

const char* pylith::meshio::DataWriterVTKDataPointsTri3::_timestepFilename = 
  "tri3_points.vtk";

const char* pylith::meshio::DataWriterVTKDataPointsTri3::_vertexFilename = 
  "tri3_points_vertex.vtk";

const PylithScalar pylith::meshio::DataWriterVTKDataPointsTri3::_time = 1.0;

const char* pylith::meshio::DataWriterVTKDataPointsTri3::_timeFormat = 
  "%3.1f";

const int pylith::meshio::DataWriterVTKDataPointsTri3::_numVertexFields = 3;
const int pylith::meshio::DataWriterVTKDataPointsTri3::_numVertices = 8;

const pylith::meshio::DataWriterData::FieldStruct
pylith::meshio::DataWriterVTKDataPointsTri3::_vertexFields[] = {
  { "displacements", topology::FieldBase::VECTOR, 2 },
  { "pressure", topology::FieldBase::SCALAR, 1 },
  { "other", topology::FieldBase::OTHER, 2 },
};
const PylithScalar pylith::meshio::DataWriterVTKDataPointsTri3::_vertexField0[] = {
  1.1, 2.2, // 0
  3.3, 4.4, // 1
  5.5, 6.6, // 2
  7.7, 8.8, // 3
  9.9, 10.0, // 4
  11.1, 12.2, // 5
  13.3, 14.4,
  15.5, 16.6,
};
const PylithScalar pylith::meshio::DataWriterVTKDataPointsTri3::_vertexField1[] = {
  2.1, // 0
  3.2, // 1
  4.3, // 2
  5.4, // 3
  6.5, // 4
  7.6, // 5
  8.7,
  9.8
};
const PylithScalar pylith::meshio::DataWriterVTKDataPointsTri3::_vertexField2[] = {
  1.2, 2.3, // 0
  3.4, 4.5, // 1
  5.6, 6.7, // 2
  7.8, 8.9, // 3
  9.0, 10.1, // 4
  11.2, 12.3, // 5
  13.4, 14.5,
  15.6, 16.7,
};

const int pylith::meshio::DataWriterVTKDataPointsTri3::_numPoints = 3;
const int pylith::meshio::DataWriterVTKDataPointsTri3::_spaceDim = 2;
const PylithScalar pylith::meshio::DataWriterVTKDataPointsTri3::_points[] = {
 -0.3333333, 0.0,
  0.0000001, 0.0,
  0.9999999, 0.0,
};


pylith::meshio::DataWriterVTKDataPointsTri3::DataWriterVTKDataPointsTri3(void)
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


pylith::meshio::DataWriterVTKDataPointsTri3::~DataWriterVTKDataPointsTri3(void)
{}


// End of file
