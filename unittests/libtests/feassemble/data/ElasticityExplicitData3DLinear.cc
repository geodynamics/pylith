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

// DO NOT EDIT THIS FILE
// This file was generated from python application elasticityexplicit.

#include "ElasticityExplicitData3DLinear.hh"

const int pylith::feassemble::ElasticityExplicitData3DLinear::_spaceDim = 3;

const int pylith::feassemble::ElasticityExplicitData3DLinear::_cellDim = 3;

const int pylith::feassemble::ElasticityExplicitData3DLinear::_numVertices = 4;

const int pylith::feassemble::ElasticityExplicitData3DLinear::_numCells = 1;

const int pylith::feassemble::ElasticityExplicitData3DLinear::_numBasis = 4;

const int pylith::feassemble::ElasticityExplicitData3DLinear::_numQuadPts = 1;

const char* pylith::feassemble::ElasticityExplicitData3DLinear::_matType = "ElasticIsotropic3D";

const char* pylith::feassemble::ElasticityExplicitData3DLinear::_matDBFilename = "data/elasticisotropic3d.spatialdb";

const int pylith::feassemble::ElasticityExplicitData3DLinear::_matId = 0;

const char* pylith::feassemble::ElasticityExplicitData3DLinear::_matLabel = "elastic isotropic 3-D";

const double pylith::feassemble::ElasticityExplicitData3DLinear::_dt =   1.00000000e-02;

const double pylith::feassemble::ElasticityExplicitData3DLinear::_vertices[] = {
 -5.00000000e-01, -1.00000000e+00, -5.00000000e-01,
  2.00000000e+00, -5.00000000e-01, -4.00000000e-01,
  1.00000000e+00, -1.00000000e-01, -3.00000000e-01,
 -2.00000000e-01,  5.00000000e-01,  2.00000000e+00,
};

const int pylith::feassemble::ElasticityExplicitData3DLinear::_cells[] = {
0,1,2,3,
};

const double pylith::feassemble::ElasticityExplicitData3DLinear::_quadPts[] = {
  2.50000000e-01,  2.50000000e-01,  2.50000000e-01,
};

const double pylith::feassemble::ElasticityExplicitData3DLinear::_quadWts[] = {
  1.66666667e-01,
};

const double pylith::feassemble::ElasticityExplicitData3DLinear::_basis[] = {
  2.50000000e-01,  2.50000000e-01,  2.50000000e-01,
  2.50000000e-01,};

const double pylith::feassemble::ElasticityExplicitData3DLinear::_basisDeriv[] = {
 -1.00000000e+00, -1.00000000e+00, -1.00000000e+00,
  1.00000000e+00,  0.00000000e+00,  0.00000000e+00,
  0.00000000e+00,  1.00000000e+00,  0.00000000e+00,
  0.00000000e+00,  0.00000000e+00,  1.00000000e+00,
};

const double pylith::feassemble::ElasticityExplicitData3DLinear::_fieldTpdt[] = {
  3.00000000e-01,  2.00000000e-01, -5.00000000e-01,
 -3.00000000e-01, -4.00000000e-01, -6.00000000e-01,
  2.00000000e-01,  6.00000000e-01,  3.00000000e-01,
 -6.00000000e-01, -1.00000000e-01, -3.00000000e-01,
};

const double pylith::feassemble::ElasticityExplicitData3DLinear::_fieldT[] = {
  8.00000000e-01,  1.00000000e-01, -6.00000000e-01,
 -1.00000000e-01, -2.00000000e-01, -5.00000000e-01,
  1.00000000e-01,  7.00000000e-01,  2.00000000e-01,
 -5.00000000e-01, -0.00000000e+00, -2.00000000e-01,
};

const double pylith::feassemble::ElasticityExplicitData3DLinear::_fieldTmdt[] = {
  1.00000000e-01,  1.00000000e-01, -3.00000000e-01,
 -2.00000000e-01, -1.00000000e-01, -5.00000000e-01,
  2.00000000e-01,  4.00000000e-01,  1.00000000e-01,
 -4.00000000e-01, -1.00000000e-01, -1.00000000e-01,
};

const double pylith::feassemble::ElasticityExplicitData3DLinear::_valsResidual[] = {
 -4.51909072e+10,  1.85615044e+10,  1.04907478e+10,
  2.74390928e+10,  8.07050438e+09,  9.68374781e+09,
  8.07109281e+09, -2.09814956e+10, -5.64925219e+09,
  9.68509281e+09, -5.64849563e+09, -1.45262522e+10,
};

const double pylith::feassemble::ElasticityExplicitData3DLinear::_valsJacobian[] = {
  8.40625000e+05,  0.00000000e+00,  0.00000000e+00,
  8.40625000e+05,  0.00000000e+00,  0.00000000e+00,
  8.40625000e+05,  0.00000000e+00,  0.00000000e+00,
  8.40625000e+05,  0.00000000e+00,  0.00000000e+00,
  0.00000000e+00,  8.40625000e+05,  0.00000000e+00,
  0.00000000e+00,  8.40625000e+05,  0.00000000e+00,
  0.00000000e+00,  8.40625000e+05,  0.00000000e+00,
  0.00000000e+00,  8.40625000e+05,  0.00000000e+00,
  0.00000000e+00,  0.00000000e+00,  8.40625000e+05,
  0.00000000e+00,  0.00000000e+00,  8.40625000e+05,
  0.00000000e+00,  0.00000000e+00,  8.40625000e+05,
  0.00000000e+00,  0.00000000e+00,  8.40625000e+05,
  8.40625000e+05,  0.00000000e+00,  0.00000000e+00,
  8.40625000e+05,  0.00000000e+00,  0.00000000e+00,
  8.40625000e+05,  0.00000000e+00,  0.00000000e+00,
  8.40625000e+05,  0.00000000e+00,  0.00000000e+00,
  0.00000000e+00,  8.40625000e+05,  0.00000000e+00,
  0.00000000e+00,  8.40625000e+05,  0.00000000e+00,
  0.00000000e+00,  8.40625000e+05,  0.00000000e+00,
  0.00000000e+00,  8.40625000e+05,  0.00000000e+00,
  0.00000000e+00,  0.00000000e+00,  8.40625000e+05,
  0.00000000e+00,  0.00000000e+00,  8.40625000e+05,
  0.00000000e+00,  0.00000000e+00,  8.40625000e+05,
  0.00000000e+00,  0.00000000e+00,  8.40625000e+05,
  8.40625000e+05,  0.00000000e+00,  0.00000000e+00,
  8.40625000e+05,  0.00000000e+00,  0.00000000e+00,
  8.40625000e+05,  0.00000000e+00,  0.00000000e+00,
  8.40625000e+05,  0.00000000e+00,  0.00000000e+00,
  0.00000000e+00,  8.40625000e+05,  0.00000000e+00,
  0.00000000e+00,  8.40625000e+05,  0.00000000e+00,
  0.00000000e+00,  8.40625000e+05,  0.00000000e+00,
  0.00000000e+00,  8.40625000e+05,  0.00000000e+00,
  0.00000000e+00,  0.00000000e+00,  8.40625000e+05,
  0.00000000e+00,  0.00000000e+00,  8.40625000e+05,
  0.00000000e+00,  0.00000000e+00,  8.40625000e+05,
  0.00000000e+00,  0.00000000e+00,  8.40625000e+05,
  8.40625000e+05,  0.00000000e+00,  0.00000000e+00,
  8.40625000e+05,  0.00000000e+00,  0.00000000e+00,
  8.40625000e+05,  0.00000000e+00,  0.00000000e+00,
  8.40625000e+05,  0.00000000e+00,  0.00000000e+00,
  0.00000000e+00,  8.40625000e+05,  0.00000000e+00,
  0.00000000e+00,  8.40625000e+05,  0.00000000e+00,
  0.00000000e+00,  8.40625000e+05,  0.00000000e+00,
  0.00000000e+00,  8.40625000e+05,  0.00000000e+00,
  0.00000000e+00,  0.00000000e+00,  8.40625000e+05,
  0.00000000e+00,  0.00000000e+00,  8.40625000e+05,
  0.00000000e+00,  0.00000000e+00,  8.40625000e+05,
  0.00000000e+00,  0.00000000e+00,  8.40625000e+05,
};

pylith::feassemble::ElasticityExplicitData3DLinear::ElasticityExplicitData3DLinear(void)
{ // constructor
  spaceDim = _spaceDim;
  cellDim = _cellDim;
  numVertices = _numVertices;
  numCells = _numCells;
  numBasis = _numBasis;
  numQuadPts = _numQuadPts;
  matType = const_cast<char*>(_matType);
  matDBFilename = const_cast<char*>(_matDBFilename);
  matId = _matId;
  matLabel = const_cast<char*>(_matLabel);
  dt = _dt;
  vertices = const_cast<double*>(_vertices);
  cells = const_cast<int*>(_cells);
  quadPts = const_cast<double*>(_quadPts);
  quadWts = const_cast<double*>(_quadWts);
  basis = const_cast<double*>(_basis);
  basisDeriv = const_cast<double*>(_basisDeriv);
  fieldTpdt = const_cast<double*>(_fieldTpdt);
  fieldT = const_cast<double*>(_fieldT);
  fieldTmdt = const_cast<double*>(_fieldTmdt);
  valsResidual = const_cast<double*>(_valsResidual);
  valsJacobian = const_cast<double*>(_valsJacobian);
} // constructor

pylith::feassemble::ElasticityExplicitData3DLinear::~ElasticityExplicitData3DLinear(void)
{}


// End of file
