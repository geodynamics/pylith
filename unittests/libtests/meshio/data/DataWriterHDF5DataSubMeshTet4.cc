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

#include "DataWriterHDF5DataSubMeshTet4.hh"

#include <assert.h> // USES assert()

const char* pylith::meshio::DataWriterHDF5DataSubMeshTet4::_meshFilename = 
  "data/tet4.mesh";

const char* pylith::meshio::DataWriterHDF5DataSubMeshTet4::_bcLabel = 
  "boundary";

const char* pylith::meshio::DataWriterHDF5DataSubMeshTet4::_faultLabel =  0;
const int pylith::meshio::DataWriterHDF5DataSubMeshTet4::_faultId = 0;

const char* pylith::meshio::DataWriterHDF5DataSubMeshTet4::_timestepFilename = 
  "tet4_surf.h5";

const char* pylith::meshio::DataWriterHDF5DataSubMeshTet4::_vertexFilename = 
  "tet4_surf_vertex.h5";

const char* pylith::meshio::DataWriterHDF5DataSubMeshTet4::_cellFilename = 
  "tet4_surf_cell.h5";

const PylithScalar pylith::meshio::DataWriterHDF5DataSubMeshTet4::_time = 1.0;

const char* pylith::meshio::DataWriterHDF5DataSubMeshTet4::_timeFormat = 
  "%3.1f";

const int pylith::meshio::DataWriterHDF5DataSubMeshTet4::_numVertices = 5;

const pylith::meshio::DataWriterData::FieldStruct
pylith::meshio::DataWriterHDF5DataSubMeshTet4::_vertexFields[4] = {
  { "pressure", topology::FieldBase::SCALAR, 1 },
  { "displacement", topology::FieldBase::VECTOR, 3 },
  { "stress", topology::FieldBase::TENSOR, 6 },
  { "other", topology::FieldBase::OTHER, 2 },
};
const PylithScalar pylith::meshio::DataWriterHDF5DataSubMeshTet4::_vertexFieldScalar[5*1] = {
  2.1, 3.2, 4.3, 5.4, 6.5,
};
const PylithScalar pylith::meshio::DataWriterHDF5DataSubMeshTet4::_vertexFieldVector[5*3] = {
  1.1, 2.2, 3.3,
  4.4, 5.5, 6.6,
  7.7, 8.8, 9.9,
  10.0, 11.1, 12.2,
  13.3, 14.4, 15.5,
};
const PylithScalar pylith::meshio::DataWriterHDF5DataSubMeshTet4::_vertexFieldTensor[5*6] = {
  1.1, 1.2, 1.3, 1.4, 1.5, 1.6,
  2.1, 2.2, 2.3, 2.4, 2.5, 2.6,
  3.1, 3.2, 3.3, 3.4, 3.5, 3.6,
  4.1, 4.2, 4.3, 4.4, 4.5, 4.6,
  5.1, 5.2, 5.3, 5.4, 5.5, 5.6,
};
const PylithScalar pylith::meshio::DataWriterHDF5DataSubMeshTet4::_vertexFieldOther[5*2] = {
  1.2, 2.3,
  3.4, 4.5,
  5.6, 6.7,
  7.8, 8.9,
  9.0, 10.1,
};

const int pylith::meshio::DataWriterHDF5DataSubMeshTet4::_numCells = 2;

const pylith::meshio::DataWriterData::FieldStruct
pylith::meshio::DataWriterHDF5DataSubMeshTet4::_cellFields[4] = {
  { "pressure", topology::FieldBase::SCALAR, 1 },
  { "traction", topology::FieldBase::VECTOR, 3 },
  { "stress", topology::FieldBase::TENSOR, 6 },
  { "other", topology::FieldBase::OTHER, 4 },
};
const PylithScalar pylith::meshio::DataWriterHDF5DataSubMeshTet4::_cellFieldScalar[2*1] = {
  2.1, 3.2
};
const PylithScalar pylith::meshio::DataWriterHDF5DataSubMeshTet4::_cellFieldVector[2*3] = {
  1.1, 2.2, 3.3,
  4.4, 5.5, 6.6,
};
const PylithScalar pylith::meshio::DataWriterHDF5DataSubMeshTet4::_cellFieldTensor[2*6] = {
  1.2, 2.3, 3.4, 4.5, 5.6, 6.7,
  7.8, 8.9, 9.0, 10.1, 11.2, 12.3,
};
const PylithScalar pylith::meshio::DataWriterHDF5DataSubMeshTet4::_cellFieldOther[2*4] = {
  1.2, 2.3, 3.4, 4.5,
  7.8, 8.9, 9.0, 10.1,
};

pylith::meshio::DataWriterHDF5DataSubMeshTet4::DataWriterHDF5DataSubMeshTet4(void)
{ // constructor
  meshFilename = const_cast<char*>(_meshFilename);
  bcLabel = const_cast<char*>(_bcLabel);
  faultLabel = const_cast<char*>(_faultLabel);
  faultId = _faultId;

  timestepFilename = const_cast<char*>(_timestepFilename);
  vertexFilename = const_cast<char*>(_vertexFilename);
  cellFilename = const_cast<char*>(_cellFilename);

  time = _time;
  timeFormat = const_cast<char*>(_timeFormat);
  
  numVertices = _numVertices;
  assert(DataWriterData::numVertexFields == numVertexFields);
  vertexFieldsInfo = const_cast<DataWriterData::FieldStruct*>(_vertexFields);
  vertexFields[0] = const_cast<PylithScalar*>(_vertexFieldScalar);
  vertexFields[1] = const_cast<PylithScalar*>(_vertexFieldVector);
  vertexFields[2] = const_cast<PylithScalar*>(_vertexFieldTensor);
  vertexFields[3] = const_cast<PylithScalar*>(_vertexFieldOther);

  numCells = _numCells;
  assert(DataWriterData::numCellFields == numCellFields);
  cellFieldsInfo = const_cast<DataWriterData::FieldStruct*>(_cellFields);
  cellFields[0] = const_cast<PylithScalar*>(_cellFieldScalar);
  cellFields[1] = const_cast<PylithScalar*>(_cellFieldVector);
  cellFields[2] = const_cast<PylithScalar*>(_cellFieldTensor);
  cellFields[3] = const_cast<PylithScalar*>(_cellFieldOther);
} // constructor

pylith::meshio::DataWriterHDF5DataSubMeshTet4::~DataWriterHDF5DataSubMeshTet4(void)
{}


// End of file
