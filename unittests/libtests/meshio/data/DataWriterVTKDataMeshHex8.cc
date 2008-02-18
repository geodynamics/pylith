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

#include "DataWriterVTKDataMeshHex8.hh"

const char* pylith::meshio::DataWriterVTKDataMeshHex8::_meshFilename = 
  "data/hex8.mesh";

const char* pylith::meshio::DataWriterVTKDataMeshHex8::_timestepFilename = 
  "hex8.vtk";

const char* pylith::meshio::DataWriterVTKDataMeshHex8::_vertexFilename = 
  "hex8_vertex.vtk";

const char* pylith::meshio::DataWriterVTKDataMeshHex8::_cellFilename = 
  "hex8_cell.vtk";

const double pylith::meshio::DataWriterVTKDataMeshHex8::_time = 1.0;

const char* pylith::meshio::DataWriterVTKDataMeshHex8::_timeFormat = 
  "%3.1f";

const int pylith::meshio::DataWriterVTKDataMeshHex8::_numVertexFields = 3;
const int pylith::meshio::DataWriterVTKDataMeshHex8::_numVertices = 12;

const pylith::meshio::DataWriterVTKData::FieldStruct
pylith::meshio::DataWriterVTKDataMeshHex8::_vertexFields[] = {
  { "displacements", VECTOR_FIELD, 3, 0 },
  { "pressure", SCALAR_FIELD, 1, 0 },
  { "other", OTHER_FIELD, 2, 0 },
};
const double pylith::meshio::DataWriterVTKDataMeshHex8::_vertexField0[] = {
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
};
const double pylith::meshio::DataWriterVTKDataMeshHex8::_vertexField1[] = {
  2.1, 3.2, 4.3, 5.4, 6.5, 7.6, 8.7, 9.8, 10.0, 12.1, 11.1, 13.1
};
const double pylith::meshio::DataWriterVTKDataMeshHex8::_vertexField2[] = {
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
  7.9, 8.0,
};

const int pylith::meshio::DataWriterVTKDataMeshHex8::_numCellFields = 3;
const int pylith::meshio::DataWriterVTKDataMeshHex8::_numCells = 1;

const pylith::meshio::DataWriterVTKData::FieldStruct
pylith::meshio::DataWriterVTKDataMeshHex8::_cellFields[] = {
  { "traction", VECTOR_FIELD, 3, 0 },
  { "pressure", SCALAR_FIELD, 1, 0 },
  { "other", TENSOR_FIELD, 6, 0 },
};
const double pylith::meshio::DataWriterVTKDataMeshHex8::_cellField0[] = {
  1.1, 2.2, 3.3,
  4.4, 5.5, 6.6
};
const double pylith::meshio::DataWriterVTKDataMeshHex8::_cellField1[] = {
  2.1, 3.2,
};
const double pylith::meshio::DataWriterVTKDataMeshHex8::_cellField2[] = {
  1.2, 2.3, 3.4, 4.5, 5.6, 6.7,
  1.1, 2.2, 3.3, 4.4, 5.5, 6.6,
};

pylith::meshio::DataWriterVTKDataMeshHex8::DataWriterVTKDataMeshHex8(void)
{ // constructor
  meshFilename = const_cast<char*>(_meshFilename);

  timestepFilename = const_cast<char*>(_timestepFilename);
  vertexFilename = const_cast<char*>(_vertexFilename);
  cellFilename = const_cast<char*>(_cellFilename);

  time = _time;
  timeFormat = const_cast<char*>(_timeFormat);
  
  numVertexFields = _numVertexFields;
  numVertices = _numVertices;
  vertexFields = const_cast<DataWriterVTKData::FieldStruct*>(_vertexFields);
  vertexFields[0].values = const_cast<double*>(_vertexField0);
  vertexFields[1].values = const_cast<double*>(_vertexField1);
  vertexFields[2].values = const_cast<double*>(_vertexField2);

  numCellFields = _numCellFields;
  numVertices = _numVertices;
  cellFields = const_cast<DataWriterVTKData::FieldStruct*>(_cellFields);
  cellFields[0].values = const_cast<double*>(_cellField0);
  cellFields[1].values = const_cast<double*>(_cellField1);
  cellFields[2].values = const_cast<double*>(_cellField2);
} // constructor

pylith::meshio::DataWriterVTKDataMeshHex8::~DataWriterVTKDataMeshHex8(void)
{}


// End of file
