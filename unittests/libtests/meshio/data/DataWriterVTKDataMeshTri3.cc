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

const char* pylith::meshio::DataWriterVTKDataMeshTri3::_meshFilename = 
  "data/tri3.mesh";

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
const int pylith::meshio::DataWriterVTKDataMeshTri3::_numVertices = 4;

const pylith::meshio::DataWriterVTKData::FieldStruct
pylith::meshio::DataWriterVTKDataMeshTri3::_vertexFields[] = {
  { "displacements", VECTOR_FIELD, 2, 0 },
  { "pressure", SCALAR_FIELD, 1, 0 },
  { "other", OTHER_FIELD, 2, 0 },
};
const double pylith::meshio::DataWriterVTKDataMeshTri3::_vertexField0[] = {
  1.1, 2.2,
  3.3, 4.4,
  5.5, 6.6,
  7.7, 8.8
};
const double pylith::meshio::DataWriterVTKDataMeshTri3::_vertexField1[] = {
  2.1, 3.2, 4.3, 5.4
};
const double pylith::meshio::DataWriterVTKDataMeshTri3::_vertexField2[] = {
  1.2, 2.3,
  3.4, 4.5,
  5.6, 6.7,
  7.8, 8.9
};

const int pylith::meshio::DataWriterVTKDataMeshTri3::_numCellFields = 3;
const int pylith::meshio::DataWriterVTKDataMeshTri3::_numCells = 2;

const pylith::meshio::DataWriterVTKData::FieldStruct
pylith::meshio::DataWriterVTKDataMeshTri3::_cellFields[] = {
  { "traction", VECTOR_FIELD, 2, 0 },
  { "pressure", SCALAR_FIELD, 1, 0 },
  { "other", TENSOR_FIELD, 3, 0 },
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

pylith::meshio::DataWriterVTKDataMeshTri3::~DataWriterVTKDataMeshTri3(void)
{}


// End of file
