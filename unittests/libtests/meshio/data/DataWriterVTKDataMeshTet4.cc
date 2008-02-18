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

#include "DataWriterVTKDataMeshTet4.hh"

const char* pylith::meshio::DataWriterVTKDataMeshTet4::_meshFilename = 
  "data/tet4.mesh";

const char* pylith::meshio::DataWriterVTKDataMeshTet4::_timestepFilename = 
  "tet4.vtk";

const char* pylith::meshio::DataWriterVTKDataMeshTet4::_vertexFilename = 
  "tet4_vertex.vtk";

const char* pylith::meshio::DataWriterVTKDataMeshTet4::_cellFilename = 
  "tet4_cell.vtk";

const double pylith::meshio::DataWriterVTKDataMeshTet4::_time = 1.0;

const char* pylith::meshio::DataWriterVTKDataMeshTet4::_timeFormat = 
  "%3.1f";

const int pylith::meshio::DataWriterVTKDataMeshTet4::_numVertexFields = 3;
const int pylith::meshio::DataWriterVTKDataMeshTet4::_numVertices = 4;

const pylith::meshio::DataWriterVTKData::FieldStruct
pylith::meshio::DataWriterVTKDataMeshTet4::_vertexFields[] = {
  { "displacements", VECTOR_FIELD, 3, 0 },
  { "pressure", SCALAR_FIELD, 1, 0 },
  { "other", OTHER_FIELD, 2, 0 },
};
const double pylith::meshio::DataWriterVTKDataMeshTet4::_vertexField0[] = {
  1.1, 2.2, 3.3,
  4.4, 5.5, 6.6,
  7.7, 8.8, 9.9,
  10.1, 11.2, 12.3
};
const double pylith::meshio::DataWriterVTKDataMeshTet4::_vertexField1[] = {
  2.1, 3.2, 4.3, 5.4,
};
const double pylith::meshio::DataWriterVTKDataMeshTet4::_vertexField2[] = {
  1.2, 2.3,
  3.4, 4.5,
  5.6, 6.7,
  7.8, 8.9,
};

const int pylith::meshio::DataWriterVTKDataMeshTet4::_numCellFields = 3;
const int pylith::meshio::DataWriterVTKDataMeshTet4::_numCells = 1;

const pylith::meshio::DataWriterVTKData::FieldStruct
pylith::meshio::DataWriterVTKDataMeshTet4::_cellFields[] = {
  { "traction", VECTOR_FIELD, 3, 0 },
  { "pressure", SCALAR_FIELD, 1, 0 },
  { "other", TENSOR_FIELD, 6, 0 },
};
const double pylith::meshio::DataWriterVTKDataMeshTet4::_cellField0[] = {
  1.1, 2.2, 3.3
};
const double pylith::meshio::DataWriterVTKDataMeshTet4::_cellField1[] = {
  2.1, 
};
const double pylith::meshio::DataWriterVTKDataMeshTet4::_cellField2[] = {
  1.2, 2.3, 3.4, 4.5, 5.6, 6.7
};

pylith::meshio::DataWriterVTKDataMeshTet4::DataWriterVTKDataMeshTet4(void)
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

pylith::meshio::DataWriterVTKDataMeshTet4::~DataWriterVTKDataMeshTet4(void)
{}


// End of file
