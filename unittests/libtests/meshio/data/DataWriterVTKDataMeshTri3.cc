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

#include "DataWriterVTKDataMeshTri3.hh"

#include <assert.h> // USES assert()

const char* pylith::meshio::DataWriterVTKDataMeshTri3::_meshFilename = 
  "data/tri3.mesh";

const char* pylith::meshio::DataWriterVTKDataMeshTri3::_faultLabel = 
  "fault";
const int pylith::meshio::DataWriterVTKDataMeshTri3::_faultId = 100;

const char* pylith::meshio::DataWriterVTKDataMeshTri3::_timestepFilename = 
  "tri3.vtk";

const char* pylith::meshio::DataWriterVTKDataMeshTri3::_vertexFilename = 
  "tri3_vertex.vtk";

const char* pylith::meshio::DataWriterVTKDataMeshTri3::_cellFilename = 
  "tri3_cell.vtk";

const double pylith::meshio::DataWriterVTKDataMeshTri3::_time = 1.0;

const char* pylith::meshio::DataWriterVTKDataMeshTri3::_timeFormat = 
  "%3.1f";

const int pylith::meshio::DataWriterVTKDataMeshTri3::_numVertexFields = 3;
const int pylith::meshio::DataWriterVTKDataMeshTri3::_numVertices = 6;

const pylith::meshio::DataWriterVTKData::FieldStruct
pylith::meshio::DataWriterVTKDataMeshTri3::_vertexFields[] = {
  { "displacements", VECTOR_FIELD, 2 },
  { "pressure", SCALAR_FIELD, 1 },
  { "other", OTHER_FIELD, 2 },
};
const double pylith::meshio::DataWriterVTKDataMeshTri3::_vertexField0[] = {
  1.1, 2.2,
  3.3, 4.4,
  5.5, 6.6,
  7.7, 8.8,
  9.9, 10.0,
  11.1, 12.2
};
const double pylith::meshio::DataWriterVTKDataMeshTri3::_vertexField1[] = {
  2.1, 3.2, 4.3, 5.4, 6.5, 7.6
};
const double pylith::meshio::DataWriterVTKDataMeshTri3::_vertexField2[] = {
  1.2, 2.3,
  3.4, 4.5,
  5.6, 6.7,
  7.8, 8.9,
  9.0, 10.1,
  11.2, 12.3
};

const int pylith::meshio::DataWriterVTKDataMeshTri3::_numCellFields = 3;
const int pylith::meshio::DataWriterVTKDataMeshTri3::_numCells = 2;

const pylith::meshio::DataWriterVTKData::FieldStruct
pylith::meshio::DataWriterVTKDataMeshTri3::_cellFields[] = {
  { "traction", VECTOR_FIELD, 2 },
  { "pressure", SCALAR_FIELD, 1 },
  { "other", TENSOR_FIELD, 3 },
};
const double pylith::meshio::DataWriterVTKDataMeshTri3::_cellField0[] = {
  1.1, 2.2,
  3.3, 4.4
};
const double pylith::meshio::DataWriterVTKDataMeshTri3::_cellField1[] = {
  2.1, 2.2,
};
const double pylith::meshio::DataWriterVTKDataMeshTri3::_cellField2[] = {
  1.2, 2.3, 3.4,
  4.5, 5.6, 6.7,
};

pylith::meshio::DataWriterVTKDataMeshTri3::DataWriterVTKDataMeshTri3(void)
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

pylith::meshio::DataWriterVTKDataMeshTri3::~DataWriterVTKDataMeshTri3(void)
{}


// End of file
