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
// Copyright (c) 2010-2014 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

#include "DataWriterHDF5DataPointsHex8.hh"

#include <assert.h> // USES assert()

const char* pylith::meshio::DataWriterHDF5DataPointsHex8::_meshFilename = 
  "data/hex8.mesh";

const char* pylith::meshio::DataWriterHDF5DataPointsHex8::_faultLabel = 
  "fault";
const int pylith::meshio::DataWriterHDF5DataPointsHex8::_faultId = 100;

const char* pylith::meshio::DataWriterHDF5DataPointsHex8::_timestepFilename = 
  "hex8_points.h5";

const char* pylith::meshio::DataWriterHDF5DataPointsHex8::_vertexFilename = 
  "hex8_points_vertex.h5";

const PylithScalar pylith::meshio::DataWriterHDF5DataPointsHex8::_time = 1.0;

const int pylith::meshio::DataWriterHDF5DataPointsHex8::_numVertices = 16;

const pylith::meshio::DataWriterData::FieldStruct
pylith::meshio::DataWriterHDF5DataPointsHex8::_vertexFields[4] = {
  { "pressure", topology::FieldBase::SCALAR, 1 },
  { "displacement", topology::FieldBase::VECTOR, 3 },
  { "stress", topology::FieldBase::TENSOR, 6 },
  { "other", topology::FieldBase::OTHER, 2 },
};
const PylithScalar pylith::meshio::DataWriterHDF5DataPointsHex8::_vertexFieldScalar[16*1] = {
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
};
const PylithScalar pylith::meshio::DataWriterHDF5DataPointsHex8::_vertexFieldVector[16*3] = {
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
};
const PylithScalar pylith::meshio::DataWriterHDF5DataPointsHex8::_vertexFieldTensor[16*6] = {
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
};
const PylithScalar pylith::meshio::DataWriterHDF5DataPointsHex8::_vertexFieldOther[16*2] = {
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
};

const int pylith::meshio::DataWriterHDF5DataPointsHex8::_numPoints = 4;
const int pylith::meshio::DataWriterHDF5DataPointsHex8::_spaceDim = 3;
const PylithScalar pylith::meshio::DataWriterHDF5DataPointsHex8::_points[4*3] = {
  -0.5, 0.0, 0.5,
  -0.00000001, 0.0, 0.0,
  -0.00000001, 0.0, 0.99999999,
   0.99999999, 0.99999999, -0.99999999,
};

const char* const pylith::meshio::DataWriterHDF5DataPointsHex8::_names[4] = {
  "ZZ.A",
  "ZZ.B",
  "ZZ.C",
  "ZZ.DD",
};


pylith::meshio::DataWriterHDF5DataPointsHex8::DataWriterHDF5DataPointsHex8(void)
{ // constructor
  meshFilename = const_cast<char*>(_meshFilename);
  faultLabel = const_cast<char*>(_faultLabel);
  faultId = _faultId;

  timestepFilename = const_cast<char*>(_timestepFilename);
  vertexFilename = const_cast<char*>(_vertexFilename);

  time = _time;
  
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
  names = _names;
} // constructor


pylith::meshio::DataWriterHDF5DataPointsHex8::~DataWriterHDF5DataPointsHex8(void)
{}


// End of file
