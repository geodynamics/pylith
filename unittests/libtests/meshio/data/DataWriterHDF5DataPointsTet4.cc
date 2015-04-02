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

#include "DataWriterHDF5DataPointsTet4.hh"

#include <assert.h> // USES assert()

const char* pylith::meshio::DataWriterHDF5DataPointsTet4::_meshFilename = 
  "data/tet4.mesh";

const char* pylith::meshio::DataWriterHDF5DataPointsTet4::_faultLabel = 
  "fault";
const int pylith::meshio::DataWriterHDF5DataPointsTet4::_faultId = 100;

const char* pylith::meshio::DataWriterHDF5DataPointsTet4::_timestepFilename = 
  "tet4_points.h5";

const char* pylith::meshio::DataWriterHDF5DataPointsTet4::_vertexFilename = 
  "tet4_points_vertex.h5";

const PylithScalar pylith::meshio::DataWriterHDF5DataPointsTet4::_time = 1.0;

const int pylith::meshio::DataWriterHDF5DataPointsTet4::_numVertices = 8;

const pylith::meshio::DataWriterData::FieldStruct
pylith::meshio::DataWriterHDF5DataPointsTet4::_vertexFields[4] = {
  { "pressure", topology::FieldBase::SCALAR, 1 },
  { "displacement", topology::FieldBase::VECTOR, 3 },
  { "stress", topology::FieldBase::TENSOR, 6 },
  { "other", topology::FieldBase::OTHER, 2 },
};
const PylithScalar pylith::meshio::DataWriterHDF5DataPointsTet4::_vertexFieldScalar[8*1] = {
  2.1, // 0
  3.2, // 1
  4.3, // 2
  5.4, // 3
  6.5, // 4
  7.6, // 5
  8.7, // 6
  9.8,  // 7
};
const PylithScalar pylith::meshio::DataWriterHDF5DataPointsTet4::_vertexFieldVector[8*3] = {
  1.1, 2.2, 3.3, // 0
  4.4, 5.5, 6.6, // 1
  7.7, 8.8, 9.9, // 2
  10.0, 11.1, 12.2, // 3
  13.3, 14.4, 15.5, // 4
  16.6, 17.7, 18.8, // 5
  19.9, 20.0, 21.1, // 6
  22.2, 23.3, 24.4, // 7
};
const PylithScalar pylith::meshio::DataWriterHDF5DataPointsTet4::_vertexFieldTensor[8*6] = {
  1.1, 1.2, 1.3, 1.4, 1.5, 1.6,
  2.1, 2.2, 2.3, 2.4, 2.5, 2.6,
  3.1, 3.2, 3.3, 3.4, 3.5, 3.6,
  4.1, 4.2, 4.3, 4.4, 4.5, 4.6,
  5.1, 5.2, 5.3, 5.4, 5.5, 5.6,
  6.1, 6.2, 6.3, 6.4, 6.5, 6.6,
  7.1, 7.2, 7.3, 7.4, 7.5, 7.6,
  8.1, 8.2, 8.3, 8.4, 8.5, 8.6,
};
const PylithScalar pylith::meshio::DataWriterHDF5DataPointsTet4::_vertexFieldOther[8*2] = {
  1.2, 2.3, // 0
  3.4, 4.5, // 1
  5.6, 6.7, // 2
  7.8, 8.9, // 3
  9.0, 10.1, // 4
  11.2, 12.3, // 5
  13.4, 14.5, // 6
  15.6, 16.7, // 7
};

const int pylith::meshio::DataWriterHDF5DataPointsTet4::_numPoints = 4;
const int pylith::meshio::DataWriterHDF5DataPointsTet4::_spaceDim = 3;
const PylithScalar pylith::meshio::DataWriterHDF5DataPointsTet4::_points[4*3] = {
  -0.33333333, 0.0, 0.33333333,
  +0.00000001, 0.0, 0.33333333,
  +0.00000001, 0.0, 0.00000001,
  0.0, -0.99999999, 0.00000001,
};

const char* const pylith::meshio::DataWriterHDF5DataPointsTet4::_names[4] = {
  "ZZ.A",
  "ZZ.B",
  "ZZ.C",
  "ZZ.DD",
};


pylith::meshio::DataWriterHDF5DataPointsTet4::DataWriterHDF5DataPointsTet4(void)
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


pylith::meshio::DataWriterHDF5DataPointsTet4::~DataWriterHDF5DataPointsTet4(void)
{}


// End of file
