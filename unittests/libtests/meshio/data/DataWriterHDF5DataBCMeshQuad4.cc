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

#include "DataWriterHDF5DataBCMeshQuad4.hh"

#include <assert.h> // USES assert()

const char* pylith::meshio::DataWriterHDF5DataBCMeshQuad4::_meshFilename = 
  "data/quad4.mesh";

const char* pylith::meshio::DataWriterHDF5DataBCMeshQuad4::_bcLabel = 
  "bc3";

const char* pylith::meshio::DataWriterHDF5DataBCMeshQuad4::_timestepFilename = 
  "quad4_bc.h5";

const char* pylith::meshio::DataWriterHDF5DataBCMeshQuad4::_vertexFilename = 
  "quad4_bc_vertex.h5";

const char* pylith::meshio::DataWriterHDF5DataBCMeshQuad4::_cellFilename = 
  "quad4_bc_cell.h5";

const PylithScalar pylith::meshio::DataWriterHDF5DataBCMeshQuad4::_time = 1.0;

const char* pylith::meshio::DataWriterHDF5DataBCMeshQuad4::_timeFormat = 
  "%3.1f";

const int pylith::meshio::DataWriterHDF5DataBCMeshQuad4::_numVertices = 3;

const pylith::meshio::DataWriterData::FieldStruct
pylith::meshio::DataWriterHDF5DataBCMeshQuad4::_vertexFields[4] = {
  { "pressure", topology::FieldBase::SCALAR, 1 },
  { "displacement", topology::FieldBase::VECTOR, 2 },
  { "stress", topology::FieldBase::TENSOR, 3 },
  { "other", topology::FieldBase::OTHER, 2 },
};
const PylithScalar pylith::meshio::DataWriterHDF5DataBCMeshQuad4::_vertexFieldScalar[3*1] = {
  2.1, 3.2, 4.3,
};
const PylithScalar pylith::meshio::DataWriterHDF5DataBCMeshQuad4::_vertexFieldVector[3*2] = {
  1.1, 2.2,
  3.3, 4.4,
  5.5, 6.6,
};
const PylithScalar pylith::meshio::DataWriterHDF5DataBCMeshQuad4::_vertexFieldTensor[3*3] = {
  1.1, 1.2, 1.3,
  2.1, 2.2, 2.3,
  3.1, 3.2, 3.3,
};
const PylithScalar pylith::meshio::DataWriterHDF5DataBCMeshQuad4::_vertexFieldOther[3*2] = {
  1.2, 2.3,
  3.4, 4.5,
  5.6, 6.7,
};

const int pylith::meshio::DataWriterHDF5DataBCMeshQuad4::_numCells = 2;

const pylith::meshio::DataWriterData::FieldStruct
pylith::meshio::DataWriterHDF5DataBCMeshQuad4::_cellFields[4] = {
  { "pressure", topology::FieldBase::SCALAR, 1 },
  { "traction", topology::FieldBase::VECTOR, 2 },
  { "stress", topology::FieldBase::TENSOR, 3 },
  { "other", topology::FieldBase::OTHER, 3 },
};
const PylithScalar pylith::meshio::DataWriterHDF5DataBCMeshQuad4::_cellFieldScalar[2*1] = {
  2.1, 3.2,
};
const PylithScalar pylith::meshio::DataWriterHDF5DataBCMeshQuad4::_cellFieldVector[2*2] = {
  1.1, 2.2,
  3.3, 4.4,
};
const PylithScalar pylith::meshio::DataWriterHDF5DataBCMeshQuad4::_cellFieldTensor[2*3] = {
  1.2, 1.3, 1.4,
  2.2, 2.3, 2.4,
};
const PylithScalar pylith::meshio::DataWriterHDF5DataBCMeshQuad4::_cellFieldOther[2*3] = {
  1.2, 2.3, 3.4,
  4.5, 5.6, 6.7,
};

pylith::meshio::DataWriterHDF5DataBCMeshQuad4::DataWriterHDF5DataBCMeshQuad4(void)
{ // constructor
  meshFilename = const_cast<char*>(_meshFilename);
  bcLabel = const_cast<char*>(_bcLabel);

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

pylith::meshio::DataWriterHDF5DataBCMeshQuad4::~DataWriterHDF5DataBCMeshQuad4(void)
{}


// End of file
