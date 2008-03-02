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

#include "DataWriterVTKDataBCMeshTri3.hh"

#include <assert.h> // USES assert()

const char* pylith::meshio::DataWriterVTKDataBCMeshTri3::_meshFilename = 
  "data/tri3.mesh";

const char* pylith::meshio::DataWriterVTKDataBCMeshTri3::_bcLabel = 
  "bc";

const char* pylith::meshio::DataWriterVTKDataBCMeshTri3::_faultLabel = 
  "fault";
const int pylith::meshio::DataWriterVTKDataBCMeshTri3::_faultId = 100;

const char* pylith::meshio::DataWriterVTKDataBCMeshTri3::_timestepFilename = 
  "tri3_bc.vtk";

const char* pylith::meshio::DataWriterVTKDataBCMeshTri3::_vertexFilename = 
  "tri3_bc_vertex.vtk";

const char* pylith::meshio::DataWriterVTKDataBCMeshTri3::_cellFilename = 
  "tri3_bc_cell.vtk";

const double pylith::meshio::DataWriterVTKDataBCMeshTri3::_time = 1.0;

const char* pylith::meshio::DataWriterVTKDataBCMeshTri3::_timeFormat = 
  "%3.1f";

const int pylith::meshio::DataWriterVTKDataBCMeshTri3::_numVertexFields = 3;
const int pylith::meshio::DataWriterVTKDataBCMeshTri3::_numVertices = 2;

const pylith::meshio::DataWriterVTKData::FieldStruct
pylith::meshio::DataWriterVTKDataBCMeshTri3::_vertexFields[] = {
  { "displacements", VECTOR_FIELD, 2 },
  { "pressure", SCALAR_FIELD, 1 },
  { "other", OTHER_FIELD, 2 },
};
const double pylith::meshio::DataWriterVTKDataBCMeshTri3::_vertexField0[] = {
  1.1, 2.2,
  3.3, 4.4,
};
const double pylith::meshio::DataWriterVTKDataBCMeshTri3::_vertexField1[] = {
  2.1, 3.2,
};
const double pylith::meshio::DataWriterVTKDataBCMeshTri3::_vertexField2[] = {
  1.2, 2.3,
  3.4, 4.5,
};

const int pylith::meshio::DataWriterVTKDataBCMeshTri3::_numCellFields = 3;
const int pylith::meshio::DataWriterVTKDataBCMeshTri3::_numCells = 1;

const pylith::meshio::DataWriterVTKData::FieldStruct
pylith::meshio::DataWriterVTKDataBCMeshTri3::_cellFields[] = {
  { "traction", VECTOR_FIELD, 2 },
  { "pressure", SCALAR_FIELD, 1 },
  { "other", TENSOR_FIELD, 3 },
};
const double pylith::meshio::DataWriterVTKDataBCMeshTri3::_cellField0[] = {
  1.1, 2.2,
};
const double pylith::meshio::DataWriterVTKDataBCMeshTri3::_cellField1[] = {
  2.1,
};
const double pylith::meshio::DataWriterVTKDataBCMeshTri3::_cellField2[] = {
  1.2, 2.3, 3.4,
};

pylith::meshio::DataWriterVTKDataBCMeshTri3::DataWriterVTKDataBCMeshTri3(void)
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

pylith::meshio::DataWriterVTKDataBCMeshTri3::~DataWriterVTKDataBCMeshTri3(void)
{}


// End of file
