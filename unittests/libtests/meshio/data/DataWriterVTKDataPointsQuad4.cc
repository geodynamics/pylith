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
// Copyright (c) 2010-2015 University of California, Davis
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

const int pylith::meshio::DataWriterVTKDataPointsQuad4::_numVertices = 6;

const pylith::meshio::DataWriterData::FieldStruct
pylith::meshio::DataWriterVTKDataPointsQuad4::_vertexFields[4] = {
  { "pressure", topology::FieldBase::SCALAR, 1 },
  { "displacement", topology::FieldBase::VECTOR, 2 },
  { "stress", topology::FieldBase::TENSOR, 3 },
  { "other", topology::FieldBase::OTHER, 2 },
};
const PylithScalar pylith::meshio::DataWriterVTKDataPointsQuad4::_vertexFieldScalar[6*1] = {
  2.1, // 0
  3.2, // 1
  4.3, // 2
  5.4, // 3
  6.5, // 4
  7.6, // 5
};
const PylithScalar pylith::meshio::DataWriterVTKDataPointsQuad4::_vertexFieldVector[6*2] = {
  1.1, 2.2, // 0
  3.3, 4.4, // 1
  5.5, 6.6, // 2
  7.7, 8.8, // 3
  9.9, 10.1, // 4
  11.2, 12.3, // 5
};
const PylithScalar pylith::meshio::DataWriterVTKDataPointsQuad4::_vertexFieldTensor[6*3] = {
  1.1, 1.2, 1.3,
  2.1, 2.2, 2.3,
  3.1, 3.2, 3.3,
  4.1, 4.2, 4.3,
  5.1, 5.2, 5.3,
  6.1, 6.2, 6.3,
};
const PylithScalar pylith::meshio::DataWriterVTKDataPointsQuad4::_vertexFieldOther[6*2] = {
  1.2, 2.3, // 0
  3.4, 4.5, // 1
  5.6, 6.7, // 2
  7.8, 8.9, // 3
  9.8, 7.6, // 4
  6.5, 5.4, // 5
};

const int pylith::meshio::DataWriterVTKDataPointsQuad4::_numPoints = 3;
const int pylith::meshio::DataWriterVTKDataPointsQuad4::_spaceDim = 2;
const PylithScalar pylith::meshio::DataWriterVTKDataPointsQuad4::_points[3*2] = {
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


pylith::meshio::DataWriterVTKDataPointsQuad4::~DataWriterVTKDataPointsQuad4(void)
{}


// End of file
