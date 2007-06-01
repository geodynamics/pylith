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
// This file was generated from python application elasticityimplicit.

#include "ElasticityImplicitData2DQuadratic.hh"

const int pylith::feassemble::ElasticityImplicitData2DQuadratic::_spaceDim = 2;

const int pylith::feassemble::ElasticityImplicitData2DQuadratic::_cellDim = 2;

const int pylith::feassemble::ElasticityImplicitData2DQuadratic::_numVertices = 6;

const int pylith::feassemble::ElasticityImplicitData2DQuadratic::_numCells = 1;

const int pylith::feassemble::ElasticityImplicitData2DQuadratic::_numBasis = 6;

const int pylith::feassemble::ElasticityImplicitData2DQuadratic::_numQuadPts = 3;

const char* pylith::feassemble::ElasticityImplicitData2DQuadratic::_matType = "ElasticPlaneStrain";

const char* pylith::feassemble::ElasticityImplicitData2DQuadratic::_matDBFilename = "data/elasticplanestrain.spatialdb";

const int pylith::feassemble::ElasticityImplicitData2DQuadratic::_matId = 0;

const char* pylith::feassemble::ElasticityImplicitData2DQuadratic::_matLabel = "elastic strain 2-D";

const double pylith::feassemble::ElasticityImplicitData2DQuadratic::_dt =   1.00000000e-02;

const double pylith::feassemble::ElasticityImplicitData2DQuadratic::_vertices[] = {
 -1.00000000e+00, -1.00000000e+00,
  1.00000000e+00,  2.00000000e-01,
 -1.50000000e+00,  5.00000000e-01,
  0.00000000e+00, -6.00000000e-01,
  2.50000000e-01,  3.50000000e-01,
 -1.25000000e+00, -2.50000000e-01,
};

const int pylith::feassemble::ElasticityImplicitData2DQuadratic::_cells[] = {
0,1,2,3,4,5,
};

const double pylith::feassemble::ElasticityImplicitData2DQuadratic::_quadPts[] = {
  6.66666667e-01,  1.66666667e-01,
  1.66666667e-01,  6.66666667e-01,
  1.66666667e-01,  1.66666667e-01,
};

const double pylith::feassemble::ElasticityImplicitData2DQuadratic::_quadWts[] = {
  1.66666667e-01,  1.66666667e-01,  1.66666667e-01,
};

const double pylith::feassemble::ElasticityImplicitData2DQuadratic::_basis[] = {
 -1.11111111e-01,  2.22222222e-01,
 -1.11111111e-01,  4.44444444e-01,
  4.44444444e-01,  1.11111111e-01,
 -1.11111111e-01, -1.11111111e-01,
  2.22222222e-01,  1.11111111e-01,
  4.44444444e-01,  4.44444444e-01,
  2.22222222e-01, -1.11111111e-01,
 -1.11111111e-01,  4.44444444e-01,
  1.11111111e-01,  4.44444444e-01,
};

const double pylith::feassemble::ElasticityImplicitData2DQuadratic::_basisDeriv[] = {
  3.33333333e-01,  3.33333333e-01,
  1.66666667e+00,  0.00000000e+00,
  0.00000000e+00, -3.33333333e-01,
 -2.00000000e+00, -2.66666667e+00,
  6.66666667e-01,  2.66666667e+00,
 -6.66666667e-01,  0.00000000e+00,
  3.33333333e-01,  3.33333333e-01,
 -3.33333333e-01,  0.00000000e+00,
  0.00000000e+00,  1.66666667e+00,
  0.00000000e+00, -6.66666667e-01,
  2.66666667e+00,  6.66666667e-01,
 -2.66666667e+00, -2.00000000e+00,
 -1.66666667e+00, -1.66666667e+00,
 -3.33333333e-01,  0.00000000e+00,
  0.00000000e+00, -3.33333333e-01,
  2.00000000e+00, -6.66666667e-01,
  6.66666667e-01,  6.66666667e-01,
 -6.66666667e-01,  2.00000000e+00,
};

const double pylith::feassemble::ElasticityImplicitData2DQuadratic::_fieldTpdt[] = {
 -4.00000000e-01, -6.00000000e-01,
  7.00000000e-01,  8.00000000e-01,
  0.00000000e+00,  2.00000000e-01,
 -5.00000000e-01, -4.00000000e-01,
  3.00000000e-01,  9.00000000e-01,
 -3.00000000e-01, -9.00000000e-01,
};

const double pylith::feassemble::ElasticityImplicitData2DQuadratic::_fieldT[] = {
 -3.00000000e-01, -4.00000000e-01,
  5.00000000e-01,  6.00000000e-01,
  0.00000000e+00,  1.00000000e-01,
 -2.00000000e-01, -3.00000000e-01,
  2.00000000e-01,  3.00000000e-01,
 -1.00000000e-01, -2.00000000e-01,
};

const double pylith::feassemble::ElasticityImplicitData2DQuadratic::_fieldTmdt[] = {
 -2.00000000e-01, -3.00000000e-01,
  3.00000000e-01,  4.00000000e-01,
  0.00000000e+00,  1.00000000e-01,
 -3.00000000e-01, -2.00000000e-01,
  1.00000000e-01,  4.00000000e-01,
 -2.00000000e-01, -6.00000000e-01,
};

const double pylith::feassemble::ElasticityImplicitData2DQuadratic::_valsResidual[] = {
 -1.48468519e+11, -1.75705556e+11,
 -2.40611111e+11, -4.21759259e+10,
 -1.05142593e+11, -3.85807407e+11,
  5.56192593e+11,  7.79718519e+11,
 -7.80303704e+11, -9.67622222e+11,
  7.18333333e+11,  7.91592593e+11,
};

const double pylith::feassemble::ElasticityImplicitData2DQuadratic::_valsJacobian[] = {
  2.09351852e+11,  8.97222222e+10,
  5.18888889e+10,  8.64814815e+09,
  1.23703704e+10,  2.47407407e+10,
 -2.31092593e+11, -3.80740741e+10,
  2.06481481e+10,  3.00000000e+09,
 -6.31666667e+10, -8.80370370e+10,
  8.97222222e+10,  2.09351852e+11,
  1.72962963e+10,  8.64814815e+09,
  1.23703704e+10,  7.42222222e+10,
 -7.68703704e+10, -3.71111111e+10,
  6.72222222e+09,  2.03703704e+09,
 -4.92407407e+10, -2.57148148e+11,
  5.18888889e+10,  1.72962963e+10,
  1.58111111e+11,  0.00000000e+00,
  0.00000000e+00, -2.29629630e+10,
 -2.10000000e+11, -6.51481481e+10,
 -3.15555556e+10,  6.51481481e+10,
  3.15555556e+10,  5.66666667e+09,
  8.64814815e+09,  8.64814815e+09,
  0.00000000e+00,  2.63518519e+10,
 -1.14814815e+10,  0.00000000e+00,
 -3.25740741e+10, -3.50000000e+10,
  3.25740741e+10, -5.25925926e+09,
  2.83333333e+09,  5.25925926e+09,
  1.23703704e+10,  1.23703704e+10,
  0.00000000e+00, -1.14814815e+10,
  4.12407407e+10,  0.00000000e+00,
 -5.90740741e+09, -8.88888889e+08,
  5.90740741e+09,  5.86296296e+10,
 -5.36111111e+10, -5.86296296e+10,
  2.47407407e+10,  7.42222222e+10,
 -2.29629630e+10,  0.00000000e+00,
  0.00000000e+00,  2.47444444e+11,
 -1.77777778e+09, -3.54444444e+10,
  1.17259259e+11,  3.54444444e+10,
 -1.17259259e+11, -3.21666667e+11,
 -2.31092593e+11, -7.68703704e+10,
 -2.10000000e+11, -3.25740741e+10,
 -5.90740741e+09, -1.77777778e+09,
  5.12629630e+11,  9.70000000e+10,
 -6.06296296e+10, -1.09296296e+11,
 -5.00000000e+09,  1.23518519e+11,
 -3.80740741e+10, -3.71111111e+10,
 -6.51481481e+10, -3.50000000e+10,
 -8.88888889e+08, -3.54444444e+10,
  9.70000000e+10,  5.01333333e+11,
 -1.21592593e+11, -4.26000000e+11,
  1.28703704e+11,  3.22222222e+10,
  2.06481481e+10,  6.72222222e+09,
 -3.15555556e+10,  3.25740741e+10,
  5.90740741e+09,  1.17259259e+11,
 -6.06296296e+10, -1.21592593e+11,
  7.22407407e+11,  1.33888889e+11,
 -6.56777778e+11, -1.68851852e+11,
  3.00000000e+09,  2.03703704e+09,
  6.51481481e+10, -5.25925926e+09,
  5.86296296e+10,  3.54444444e+10,
 -1.09296296e+11, -4.26000000e+11,
  1.33888889e+11,  5.36296296e+11,
 -1.51370370e+11, -1.42518519e+11,
 -6.31666667e+10, -4.92407407e+10,
  3.15555556e+10,  2.83333333e+09,
 -5.36111111e+10, -1.17259259e+11,
 -5.00000000e+09,  1.28703704e+11,
 -6.56777778e+11, -1.51370370e+11,
  7.47000000e+11,  1.86333333e+11,
 -8.80370370e+10, -2.57148148e+11,
  5.66666667e+09,  5.25925926e+09,
 -5.86296296e+10, -3.21666667e+11,
  1.23518519e+11,  3.22222222e+10,
 -1.68851852e+11, -1.42518519e+11,
  1.86333333e+11,  6.83851852e+11,
};

pylith::feassemble::ElasticityImplicitData2DQuadratic::ElasticityImplicitData2DQuadratic(void)
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

pylith::feassemble::ElasticityImplicitData2DQuadratic::~ElasticityImplicitData2DQuadratic(void)
{}


// End of file
