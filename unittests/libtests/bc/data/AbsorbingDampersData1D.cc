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

#include "AbsorbingDampersData1D.hh"

const char* pylith::bc::AbsorbingDampersData1D::_meshFilename = "data/meshBC1D.txt";

const int pylith::bc::AbsorbingDampersData1D::_spaceDim = 1;

const int pylith::bc::AbsorbingDampersData1D::_cellDim = 0;

const int pylith::bc::AbsorbingDampersData1D::_numVertices = 2;

const int pylith::bc::AbsorbingDampersData1D::_numCells = 1;

const int pylith::bc::AbsorbingDampersData1D::_numBasis = 2;

const int pylith::bc::AbsorbingDampersData1D::_numQuadPts = 1;

const char* pylith::bc::AbsorbingDampersData1D::_spatialDBFilename = "data/elasticstrain1d.spatialdb";

const double pylith::bc::AbsorbingDampersData1D::_dt =   1.00000000e-02;

const double pylith::bc::AbsorbingDampersData1D::_vertices[] = {
 -2.50000000e-01,
  2.00000000e+00,
};

const int pylith::bc::AbsorbingDampersData1D::_cells[] = {
0,1,
};

const double pylith::bc::AbsorbingDampersData1D::_verticesRef[] = {
 -1.00000000e+00,
  1.00000000e+00,
};

const double pylith::bc::AbsorbingDampersData1D::_quadPts[] = {
  0.00000000e+00,
};

const double pylith::bc::AbsorbingDampersData1D::_quadWts[] = {
  2.00000000e+00,
};

const double pylith::bc::AbsorbingDampersData1D::_basis[] = {
  5.00000000e-01,
  5.00000000e-01,
};

const double pylith::bc::AbsorbingDampersData1D::_basisDerivRef[] = {
 -5.00000000e-01,
  5.00000000e-01,
};

const double pylith::bc::AbsorbingDampersData1D::_fieldTpdt[] = {
  1.20000000e+00,
  1.70000000e+00,
};

const double pylith::bc::AbsorbingDampersData1D::_fieldT[] = {
  1.10000000e+00,
  1.50000000e+00,
};

const double pylith::bc::AbsorbingDampersData1D::_fieldTmdt[] = {
  1.00000000e+00,
  1.30000000e+00,
};

const double pylith::bc::AbsorbingDampersData1D::_valsResidual[] = {
  1.60000000e+10,
 -1.60000000e+10,
};

const double pylith::bc::AbsorbingDampersData1D::_valsJacobian[] = {
  1.40625000e+07,
  1.40625000e+07,
  1.40625000e+07,
  1.40625000e+07,
};

pylith::bc::AbsorbingDampersData1D::AbsorbingDampersData1D(void)
{ // constructor
  meshFilename = const_cast<char*>(_meshFilename);
  spaceDim = _spaceDim;
  cellDim = _cellDim;
  numVertices = _numVertices;
  numCells = _numCells;
  numBasis = _numBasis;
  numQuadPts = _numQuadPts;
  spatialDBFilename = const_cast<char*>(_spatialDBFilename);
  dt = _dt;
  vertices = const_cast<double*>(_vertices);
  cells = const_cast<int*>(_cells);
  verticesRef = const_cast<double*>(_verticesRef);
  quadPts = const_cast<double*>(_quadPts);
  quadWts = const_cast<double*>(_quadWts);
  basis = const_cast<double*>(_basis);
  basisDerivRef = const_cast<double*>(_basisDerivRef);
  fieldTpdt = const_cast<double*>(_fieldTpdt);
  fieldT = const_cast<double*>(_fieldT);
  fieldTmdt = const_cast<double*>(_fieldTmdt);
  valsResidual = const_cast<double*>(_valsResidual);
  valsJacobian = const_cast<double*>(_valsJacobian);
} // constructor

pylith::bc::AbsorbingDampersData1D::~AbsorbingDampersData1D(void)
{}


// End of file
