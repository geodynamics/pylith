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

#include "DataWriterVTKDataMeshLine2.hh"

#include <assert.h> // USES assert()

const char* pylith::meshio::DataWriterVTKDataMeshLine2::_meshFilename = 
  "data/line2.mesh";

const char* pylith::meshio::DataWriterVTKDataMeshLine2::_timestepFilename = 
  "line2.vtk";

const char* pylith::meshio::DataWriterVTKDataMeshLine2::_vertexFilename = 
  "line2_vertex.vtk";

const char* pylith::meshio::DataWriterVTKDataMeshLine2::_cellFilename = 
  "line2_cell.vtk";

const double pylith::meshio::DataWriterVTKDataMeshLine2::_time = 1.0;

const char* pylith::meshio::DataWriterVTKDataMeshLine2::_timeFormat = 
  "%3.1f";

const int pylith::meshio::DataWriterVTKDataMeshLine2::_numVertexFields = 3;
const int pylith::meshio::DataWriterVTKDataMeshLine2::_numVertices = 3;

const pylith::meshio::DataWriterVTKData::FieldStruct
pylith::meshio::DataWriterVTKDataMeshLine2::_vertexFields[] = {
  { "displacements", VECTOR_FIELD, 1 },
  { "pressure", SCALAR_FIELD, 1 },
  { "other", OTHER_FIELD, 2 },
};
const double pylith::meshio::DataWriterVTKDataMeshLine2::_vertexField0[] = {
  1.1, 2.2, 3.3,
};
const double pylith::meshio::DataWriterVTKDataMeshLine2::_vertexField1[] = {
  2.1, 3.2, 4.3,
};
const double pylith::meshio::DataWriterVTKDataMeshLine2::_vertexField2[] = {
  1.2, 2.3,
  3.4, 4.5,
  5.6, 6.7,
};

const int pylith::meshio::DataWriterVTKDataMeshLine2::_numCellFields = 3;
const int pylith::meshio::DataWriterVTKDataMeshLine2::_numCells = 2;

const pylith::meshio::DataWriterVTKData::FieldStruct
pylith::meshio::DataWriterVTKDataMeshLine2::_cellFields[] = {
  { "traction", VECTOR_FIELD, 1 },
  { "pressure", SCALAR_FIELD, 1 },
  { "other", TENSOR_FIELD, 1 },
};
const double pylith::meshio::DataWriterVTKDataMeshLine2::_cellField0[] = {
  1.1, 2.2,
};
const double pylith::meshio::DataWriterVTKDataMeshLine2::_cellField1[] = {
  2.1, 2.2,
};
const double pylith::meshio::DataWriterVTKDataMeshLine2::_cellField2[] = {
  1.2, 2.3,
};

pylith::meshio::DataWriterVTKDataMeshLine2::DataWriterVTKDataMeshLine2(void)
{ // constructor
  meshFilename = const_cast<char*>(_meshFilename);

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
  numVertices = _numVertices;
  assert(3 == numCellFields);
  cellFieldsInfo = const_cast<DataWriterVTKData::FieldStruct*>(_cellFields);
  cellFields[0] = const_cast<double*>(_cellField0);
  cellFields[1] = const_cast<double*>(_cellField1);
  cellFields[2] = const_cast<double*>(_cellField2);
} // constructor

pylith::meshio::DataWriterVTKDataMeshLine2::~DataWriterVTKDataMeshLine2(void)
{}


// End of file
