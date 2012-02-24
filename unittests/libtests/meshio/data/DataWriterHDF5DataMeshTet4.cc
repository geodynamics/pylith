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

#include "DataWriterHDF5DataMeshTet4.hh"

#include <assert.h> // USES assert()

const char* pylith::meshio::DataWriterHDF5DataMeshTet4::_meshFilename = 
  "data/tet4.mesh";

const char* pylith::meshio::DataWriterHDF5DataMeshTet4::_faultLabel = 
  "fault";
const int pylith::meshio::DataWriterHDF5DataMeshTet4::_faultId = 100;

const char* pylith::meshio::DataWriterHDF5DataMeshTet4::_timestepFilename = 
  "tet4.h5";

const char* pylith::meshio::DataWriterHDF5DataMeshTet4::_vertexFilename = 
  "tet4_vertex.h5";

const char* pylith::meshio::DataWriterHDF5DataMeshTet4::_cellFilename = 
  "tet4_cell.h5";

const double pylith::meshio::DataWriterHDF5DataMeshTet4::_time = 1.0;

const char* pylith::meshio::DataWriterHDF5DataMeshTet4::_timeFormat = 
  "%3.1f";

const int pylith::meshio::DataWriterHDF5DataMeshTet4::_numVertexFields = 3;
const int pylith::meshio::DataWriterHDF5DataMeshTet4::_numVertices = 11;

const pylith::meshio::DataWriterData::FieldStruct
pylith::meshio::DataWriterHDF5DataMeshTet4::_vertexFields[] = {
  { "displacements", topology::FieldBase::VECTOR, 3 },
  { "pressure", topology::FieldBase::SCALAR, 1 },
  { "other", topology::FieldBase::OTHER, 2 },
};
const double pylith::meshio::DataWriterHDF5DataMeshTet4::_vertexField0[] = {
  1.1, 2.2, 3.3,
  4.4, 5.5, 6.6,
  7.7, 8.8, 9.9,
  10.0, 11.1, 12.2,
  13.3, 14.4, 15.5,
  16.6, 17.7, 18.8,
  19.9, 20.0, 21.1,
  22.2, 23.3, 24.4,
  25.5, 26.6, 27.7,
  28.8, 29.9, 30.0,
  31.1, 32.2, 33.3,
};
const double pylith::meshio::DataWriterHDF5DataMeshTet4::_vertexField1[] = {
  2.1, 3.2, 4.3, 5.4, 6.5, 7.6, 8.7, 9.8, 10.9, 11.0, 12.1,
};
const double pylith::meshio::DataWriterHDF5DataMeshTet4::_vertexField2[] = {
  1.2, 2.3,
  3.4, 4.5,
  5.6, 6.7,
  7.8, 8.9,
  9.0, 10.1,
  11.2, 12.3,
  13.4, 14.5,
  15.6, 16.7,
  17.8, 18.9,
  19.0, 20.1,
  21.2, 22.3,
};

const int pylith::meshio::DataWriterHDF5DataMeshTet4::_numCellFields = 3;
const int pylith::meshio::DataWriterHDF5DataMeshTet4::_numCells = 3;

const pylith::meshio::DataWriterData::FieldStruct
pylith::meshio::DataWriterHDF5DataMeshTet4::_cellFields[] = {
  { "traction", topology::FieldBase::VECTOR, 3 },
  { "pressure", topology::FieldBase::SCALAR, 1 },
  { "other", topology::FieldBase::TENSOR, 6 },
};
const double pylith::meshio::DataWriterHDF5DataMeshTet4::_cellField0[] = {
  1.1, 2.2, 3.3,
  4.4, 5.5, 6.6,
  7.7, 8.8, 9.9
};
const double pylith::meshio::DataWriterHDF5DataMeshTet4::_cellField1[] = {
  2.1, 3.2, 4.3
};
const double pylith::meshio::DataWriterHDF5DataMeshTet4::_cellField2[] = {
  1.2, 2.3, 3.4, 4.5, 5.6, 6.7,
  7.8, 8.9, 9.0, 10.1, 11.2, 12.3,
  13.4, 14.5, 15.6, 16.7, 17.8, 18.9
};

pylith::meshio::DataWriterHDF5DataMeshTet4::DataWriterHDF5DataMeshTet4(void)
{ // constructor
  meshFilename = const_cast<char*>(_meshFilename);
  faultLabel = const_cast<char*>(_faultLabel);
  faultId = _faultId;

  timestepFilename = const_cast<char*>(_timestepFilename);
  vertexFilename = const_cast<char*>(_vertexFilename);
  cellFilename = const_cast<char*>(_cellFilename);

  time = _time;
  timeFormat = const_cast<char*>(_timeFormat);
  
  numVertexFields = _numVertexFields;
  assert(3 == numVertexFields);
  numVertices = _numVertices;
  vertexFieldsInfo = const_cast<DataWriterData::FieldStruct*>(_vertexFields);
  vertexFields[0] = const_cast<double*>(_vertexField0);
  vertexFields[1] = const_cast<double*>(_vertexField1);
  vertexFields[2] = const_cast<double*>(_vertexField2);

  numCellFields = _numCellFields;
  assert(3 == numCellFields);
  numCells = _numCells;
  cellFieldsInfo = const_cast<DataWriterData::FieldStruct*>(_cellFields);
  cellFields[0] = const_cast<double*>(_cellField0);
  cellFields[1] = const_cast<double*>(_cellField1);
  cellFields[2] = const_cast<double*>(_cellField2);
} // constructor

pylith::meshio::DataWriterHDF5DataMeshTet4::~DataWriterHDF5DataMeshTet4(void)
{}


// End of file
