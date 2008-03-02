// -*- C++ -*-
//
// ======================================================================
//
//                           Brad T. Aagaard
//                        U.S. Geological Survey
//
// {LicenseText}
//
// ======================================================================
//

#include "DataWriterVTKDataSubMeshHex8.hh"

#include <assert.h> // USES assert()

const char* pylith::meshio::DataWriterVTKDataSubMeshHex8::_meshFilename = 
  "data/hex8.mesh";

const char* pylith::meshio::DataWriterVTKDataSubMeshHex8::_cellsLabel = 
  "material-id";
const int pylith::meshio::DataWriterVTKDataSubMeshHex8::_labelId = 0;

const char* pylith::meshio::DataWriterVTKDataSubMeshHex8::_faultLabel = 
  "fault";
const int pylith::meshio::DataWriterVTKDataSubMeshHex8::_faultId = 100;

const char* pylith::meshio::DataWriterVTKDataSubMeshHex8::_timestepFilename = 
  "hex8_sub.vtk";

const char* pylith::meshio::DataWriterVTKDataSubMeshHex8::_vertexFilename = 
  "hex8_sub_vertex.vtk";

const char* pylith::meshio::DataWriterVTKDataSubMeshHex8::_cellFilename = 
  "hex8_sub_cell.vtk";

const double pylith::meshio::DataWriterVTKDataSubMeshHex8::_time = 1.0;

const char* pylith::meshio::DataWriterVTKDataSubMeshHex8::_timeFormat = 
  "%3.1f";

const int pylith::meshio::DataWriterVTKDataSubMeshHex8::_numVertexFields = 3;
const int pylith::meshio::DataWriterVTKDataSubMeshHex8::_numVertices = 20;

const pylith::meshio::DataWriterVTKData::FieldStruct
pylith::meshio::DataWriterVTKDataSubMeshHex8::_vertexFields[] = {
  { "displacements", VECTOR_FIELD, 3 },
  { "pressure", SCALAR_FIELD, 1 },
  { "other", OTHER_FIELD, 2 },
};
const double pylith::meshio::DataWriterVTKDataSubMeshHex8::_vertexField0[] = {
  1.1, 2.2, 3.3,
  4.4, 5.5, 6.6,
  7.7, 8.8, 9.9,
  10.1, 11.2, 12.3,
  1.2, 2.3, 3.4,
  4.5, 5.6, 6.7,
  7.8, 8.9, 9.0,
  10.2, 11.3, 12.4,
  1.3, 2.4, 3.5,
  4.6, 5.7, 6.8,
  7.9, 8.0, 9.1,
  10.2, 11.3, 12.4,
  13.5, 14.6, 15.7,
  16.8, 17.9, 18.1,
  19.2, 20.3, 21.4,
  22.5, 23.6, 24.7,
  25.5, 26.6, 27.7,
  28.8, 29.9, 30.0,
  31.1, 32.2, 33.3,
  34.4, 35.5, 36.6
};
const double pylith::meshio::DataWriterVTKDataSubMeshHex8::_vertexField1[] = {
  2.1, 3.2, 4.3, 5.4, 6.5, 7.6, 8.7, 9.8,
  10.0, 12.1, 11.1, 13.1, 14.1, 15.1, 16.1, 17.1,
  18.1, 19.1, 20.1, 21.1
};
const double pylith::meshio::DataWriterVTKDataSubMeshHex8::_vertexField2[] = {
  1.2, 2.3,
  3.4, 4.5,
  5.6, 6.7,
  7.8, 8.9,
  1.3, 2.4,
  3.5, 4.6,
  5.7, 6.8,
  7.9, 8.0,
  1.3, 2.4,
  3.5, 4.6,
  5.7, 6.8,
  8.0, 1.4,
  2.5, 3.6,
  4.8, 1.5,
  2.6, 3.7,
  4.8, 5.9,
  1.1, 1.2,
  2.1, 2.2,
  3.1, 3.2,
  4.1, 4.2,
};

const int pylith::meshio::DataWriterVTKDataSubMeshHex8::_numCellFields = 3;
const int pylith::meshio::DataWriterVTKDataSubMeshHex8::_numCells = 1;

const pylith::meshio::DataWriterVTKData::FieldStruct
pylith::meshio::DataWriterVTKDataSubMeshHex8::_cellFields[] = {
  { "traction", VECTOR_FIELD, 3 },
  { "pressure", SCALAR_FIELD, 1 },
  { "other", TENSOR_FIELD, 6 },
};
const double pylith::meshio::DataWriterVTKDataSubMeshHex8::_cellField0[] = {
  1.1, 2.2, 3.3,
};
const double pylith::meshio::DataWriterVTKDataSubMeshHex8::_cellField1[] = {
  2.1,
};
const double pylith::meshio::DataWriterVTKDataSubMeshHex8::_cellField2[] = {
  1.2, 2.3, 3.4, 4.5, 5.6, 6.7,
};

pylith::meshio::DataWriterVTKDataSubMeshHex8::DataWriterVTKDataSubMeshHex8(void)
{ // constructor
  meshFilename = const_cast<char*>(_meshFilename);
  cellsLabel = const_cast<char*>(_cellsLabel);
  labelId = _labelId;
  faultLabel = const_cast<char*>(_faultLabel);
  faultId = _faultId;

  timestepFilename = const_cast<char*>(_timestepFilename);
  vertexFilename = const_cast<char*>(_vertexFilename);
  cellFilename = const_cast<char*>(_cellFilename);

  time = _time;
  timeFormat = const_cast<char*>(_timeFormat);
  
  numVertexFields = _numVertexFields;
  numVertices = _numVertices;
  assert(3 == numVertexFields);
  vertexFieldsInfo = const_cast<DataWriterVTKData::FieldStruct*>(_vertexFields);
  vertexFields[0] = const_cast<double*>(_vertexField0);
  vertexFields[1] = const_cast<double*>(_vertexField1);
  vertexFields[2] = const_cast<double*>(_vertexField2);

  numCellFields = _numCellFields;
  numCells = _numCells;
  assert(3 == numCellFields);
  cellFieldsInfo = const_cast<DataWriterVTKData::FieldStruct*>(_cellFields);
  cellFields[0] = const_cast<double*>(_cellField0);
  cellFields[1] = const_cast<double*>(_cellField1);
  cellFields[2] = const_cast<double*>(_cellField2);
} // constructor

pylith::meshio::DataWriterVTKDataSubMeshHex8::~DataWriterVTKDataSubMeshHex8(void)
{}


// End of file
