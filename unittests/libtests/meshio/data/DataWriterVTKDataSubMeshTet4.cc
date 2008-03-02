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

#include "DataWriterVTKDataSubMeshTet4.hh"

#include <assert.h> // USES assert()

const char* pylith::meshio::DataWriterVTKDataSubMeshTet4::_meshFilename = 
  "data/tet4.mesh";

const char* pylith::meshio::DataWriterVTKDataSubMeshTet4::_cellsLabel = 
  "material-id";
const int pylith::meshio::DataWriterVTKDataSubMeshTet4::_labelId = 1;

const char* pylith::meshio::DataWriterVTKDataSubMeshTet4::_faultLabel = 
  "fault";
const int pylith::meshio::DataWriterVTKDataSubMeshTet4::_faultId = 100;

const char* pylith::meshio::DataWriterVTKDataSubMeshTet4::_timestepFilename = 
  "tet4_sub.vtk";

const char* pylith::meshio::DataWriterVTKDataSubMeshTet4::_vertexFilename = 
  "tet4_sub_vertex.vtk";

const char* pylith::meshio::DataWriterVTKDataSubMeshTet4::_cellFilename = 
  "tet4_sub_cell.vtk";

const double pylith::meshio::DataWriterVTKDataSubMeshTet4::_time = 1.0;

const char* pylith::meshio::DataWriterVTKDataSubMeshTet4::_timeFormat = 
  "%3.1f";

const int pylith::meshio::DataWriterVTKDataSubMeshTet4::_numVertexFields = 3;
const int pylith::meshio::DataWriterVTKDataSubMeshTet4::_numVertices = 11;

const pylith::meshio::DataWriterVTKData::FieldStruct
pylith::meshio::DataWriterVTKDataSubMeshTet4::_vertexFields[] = {
  { "displacements", VECTOR_FIELD, 3 },
  { "pressure", SCALAR_FIELD, 1 },
  { "other", OTHER_FIELD, 2 },
};
const double pylith::meshio::DataWriterVTKDataSubMeshTet4::_vertexField0[] = {
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
const double pylith::meshio::DataWriterVTKDataSubMeshTet4::_vertexField1[] = {
  2.1, 3.2, 4.3, 5.4, 6.5, 7.6, 8.7, 9.8, 10.9, 11.0, 12.1
};
const double pylith::meshio::DataWriterVTKDataSubMeshTet4::_vertexField2[] = {
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
  21.2, 22.3
};

const int pylith::meshio::DataWriterVTKDataSubMeshTet4::_numCellFields = 3;
const int pylith::meshio::DataWriterVTKDataSubMeshTet4::_numCells = 2;

const pylith::meshio::DataWriterVTKData::FieldStruct
pylith::meshio::DataWriterVTKDataSubMeshTet4::_cellFields[] = {
  { "traction", VECTOR_FIELD, 3 },
  { "pressure", SCALAR_FIELD, 1 },
  { "other", TENSOR_FIELD, 6 },
};
const double pylith::meshio::DataWriterVTKDataSubMeshTet4::_cellField0[] = {
  1.1, 2.2, 3.3,
  4.4, 5.5, 6.6,
};
const double pylith::meshio::DataWriterVTKDataSubMeshTet4::_cellField1[] = {
  2.1, 3.2
};
const double pylith::meshio::DataWriterVTKDataSubMeshTet4::_cellField2[] = {
  1.2, 2.3, 3.4, 4.5, 5.6, 6.7,
  7.8, 8.9, 9.0, 10.1, 11.2, 12.3,
};

pylith::meshio::DataWriterVTKDataSubMeshTet4::DataWriterVTKDataSubMeshTet4(void)
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
  assert(3 == numVertexFields);
  numVertices = _numVertices;
  vertexFieldsInfo = const_cast<DataWriterVTKData::FieldStruct*>(_vertexFields);
  vertexFields[0] = const_cast<double*>(_vertexField0);
  vertexFields[1] = const_cast<double*>(_vertexField1);
  vertexFields[2] = const_cast<double*>(_vertexField2);

  numCellFields = _numCellFields;
  assert(3 == numCellFields);
  numCells = _numCells;
  cellFieldsInfo = const_cast<DataWriterVTKData::FieldStruct*>(_cellFields);
  cellFields[0] = const_cast<double*>(_cellField0);
  cellFields[1] = const_cast<double*>(_cellField1);
  cellFields[2] = const_cast<double*>(_cellField2);
} // constructor

pylith::meshio::DataWriterVTKDataSubMeshTet4::~DataWriterVTKDataSubMeshTet4(void)
{}


// End of file
