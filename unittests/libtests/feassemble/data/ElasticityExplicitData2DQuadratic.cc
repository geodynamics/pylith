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

#include "ElasticityExplicitData2DQuadratic.hh"

const int pylith::feassemble::ElasticityExplicitData2DQuadratic::_spaceDim = 2;

const int pylith::feassemble::ElasticityExplicitData2DQuadratic::_cellDim = 2;

const int pylith::feassemble::ElasticityExplicitData2DQuadratic::_numVertices = 6;

const int pylith::feassemble::ElasticityExplicitData2DQuadratic::_numCells = 1;

const int pylith::feassemble::ElasticityExplicitData2DQuadratic::_numBasis = 6;

const int pylith::feassemble::ElasticityExplicitData2DQuadratic::_numQuadPts = 3;

const char* pylith::feassemble::ElasticityExplicitData2DQuadratic::_matType = "ElasticPlaneStrain";

const char* pylith::feassemble::ElasticityExplicitData2DQuadratic::_matDBFilename = "data/elasticplanestrain.spatialdb";

const int pylith::feassemble::ElasticityExplicitData2DQuadratic::_matId = 0;

const char* pylith::feassemble::ElasticityExplicitData2DQuadratic::_matLabel = "elastic strain 2-D";

const double pylith::feassemble::ElasticityExplicitData2DQuadratic::_dt =   1.00000000e-02;

const double pylith::feassemble::ElasticityExplicitData2DQuadratic::_vertices[] = {
 -1.00000000e+00, -1.00000000e+00,
  1.00000000e+00,  2.00000000e-01,
 -1.50000000e+00,  5.00000000e-01,
  0.00000000e+00, -6.00000000e-01,
  2.50000000e-01,  3.50000000e-01,
 -1.25000000e+00, -2.50000000e-01,
};

const int pylith::feassemble::ElasticityExplicitData2DQuadratic::_cells[] = {
0,1,2,3,4,5,
};

const double pylith::feassemble::ElasticityExplicitData2DQuadratic::_quadPts[] = {
  6.66666667e-01,  1.66666667e-01,
  1.66666667e-01,  6.66666667e-01,
  1.66666667e-01,  1.66666667e-01,
};

const double pylith::feassemble::ElasticityExplicitData2DQuadratic::_quadWts[] = {
  1.66666667e-01,  1.66666667e-01,  1.66666667e-01,
};

const double pylith::feassemble::ElasticityExplicitData2DQuadratic::_basis[] = {
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

const double pylith::feassemble::ElasticityExplicitData2DQuadratic::_basisDeriv[] = {
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

const double pylith::feassemble::ElasticityExplicitData2DQuadratic::_fieldTpdt[] = {
 -4.00000000e-01, -6.00000000e-01,
  7.00000000e-01,  8.00000000e-01,
  0.00000000e+00,  2.00000000e-01,
 -5.00000000e-01, -4.00000000e-01,
  3.00000000e-01,  9.00000000e-01,
 -3.00000000e-01, -9.00000000e-01,
};

const double pylith::feassemble::ElasticityExplicitData2DQuadratic::_fieldT[] = {
 -3.00000000e-01, -4.00000000e-01,
  5.00000000e-01,  6.00000000e-01,
  0.00000000e+00,  1.00000000e-01,
 -2.00000000e-01, -3.00000000e-01,
  2.00000000e-01,  3.00000000e-01,
 -1.00000000e-01, -2.00000000e-01,
};

const double pylith::feassemble::ElasticityExplicitData2DQuadratic::_fieldTmdt[] = {
 -2.00000000e-01, -3.00000000e-01,
  3.00000000e-01,  4.00000000e-01,
  0.00000000e+00,  1.00000000e-01,
 -3.00000000e-01, -2.00000000e-01,
  1.00000000e-01,  4.00000000e-01,
 -2.00000000e-01, -6.00000000e-01,
};

const double pylith::feassemble::ElasticityExplicitData2DQuadratic::_valsResidual[] = {
 -2.15070691e+10,  4.58302555e+09,
 -1.24763329e+11, -3.03758597e+10,
 -2.17573798e+10, -9.72849436e+10,
  2.25394811e+11,  3.72087616e+11,
 -3.39305499e+11, -4.12213687e+11,
  2.81950760e+11,  1.63211865e+11,
};

const double pylith::feassemble::ElasticityExplicitData2DQuadratic::_valsJacobian[] = {
  1.27743484e+06,  0.00000000e+00,
 -4.66392318e+05,  0.00000000e+00,
 -8.11042524e+05,  0.00000000e+00,
  6.31001372e+05,  0.00000000e+00,
 -1.45747599e+06,  0.00000000e+00,
  2.86351166e+05,  0.00000000e+00,
  0.00000000e+00,  1.27743484e+06,
  0.00000000e+00, -4.66392318e+05,
  0.00000000e+00, -8.11042524e+05,
  0.00000000e+00,  6.31001372e+05,
  0.00000000e+00, -1.45747599e+06,
  0.00000000e+00,  2.86351166e+05,
 -4.66392318e+05,  0.00000000e+00,
  1.19513032e+06,  0.00000000e+00,
 -7.28737997e+05,  0.00000000e+00,
  3.01783265e+05,  0.00000000e+00,
  3.94375857e+04,  0.00000000e+00,
 -1.62208505e+06,  0.00000000e+00,
  0.00000000e+00, -4.66392318e+05,
  0.00000000e+00,  1.19513032e+06,
  0.00000000e+00, -7.28737997e+05,
  0.00000000e+00,  3.01783265e+05,
  0.00000000e+00,  3.94375857e+04,
  0.00000000e+00, -1.62208505e+06,
 -8.11042524e+05,  0.00000000e+00,
 -7.28737997e+05,  0.00000000e+00,
  1.53978052e+06,  0.00000000e+00,
 -9.32784636e+05,  0.00000000e+00,
  1.41803841e+06,  0.00000000e+00,
  1.33573388e+06,  0.00000000e+00,
  0.00000000e+00, -8.11042524e+05,
  0.00000000e+00, -7.28737997e+05,
  0.00000000e+00,  1.53978052e+06,
  0.00000000e+00, -9.32784636e+05,
  0.00000000e+00,  1.41803841e+06,
  0.00000000e+00,  1.33573388e+06,
  6.31001372e+05,  0.00000000e+00,
  3.01783265e+05,  0.00000000e+00,
 -9.32784636e+05,  0.00000000e+00,
  6.34430727e+06,  0.00000000e+00,
  4.78052126e+06,  0.00000000e+00,
  5.10973937e+06,  0.00000000e+00,
  0.00000000e+00,  6.31001372e+05,
  0.00000000e+00,  3.01783265e+05,
  0.00000000e+00, -9.32784636e+05,
  0.00000000e+00,  6.34430727e+06,
  0.00000000e+00,  4.78052126e+06,
  0.00000000e+00,  5.10973937e+06,
 -1.45747599e+06,  0.00000000e+00,
  3.94375857e+04,  0.00000000e+00,
  1.41803841e+06,  0.00000000e+00,
  4.78052126e+06,  0.00000000e+00,
  7.65603567e+06,  0.00000000e+00,
  6.15912209e+06,  0.00000000e+00,
  0.00000000e+00, -1.45747599e+06,
  0.00000000e+00,  3.94375857e+04,
  0.00000000e+00,  1.41803841e+06,
  0.00000000e+00,  4.78052126e+06,
  0.00000000e+00,  7.65603567e+06,
  0.00000000e+00,  6.15912209e+06,
  2.86351166e+05,  0.00000000e+00,
 -1.62208505e+06,  0.00000000e+00,
  1.33573388e+06,  0.00000000e+00,
  5.10973937e+06,  0.00000000e+00,
  6.15912209e+06,  0.00000000e+00,
  8.06755830e+06,  0.00000000e+00,
  0.00000000e+00,  2.86351166e+05,
  0.00000000e+00, -1.62208505e+06,
  0.00000000e+00,  1.33573388e+06,
  0.00000000e+00,  5.10973937e+06,
  0.00000000e+00,  6.15912209e+06,
  0.00000000e+00,  8.06755830e+06,
};

pylith::feassemble::ElasticityExplicitData2DQuadratic::ElasticityExplicitData2DQuadratic(void)
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

pylith::feassemble::ElasticityExplicitData2DQuadratic::~ElasticityExplicitData2DQuadratic(void)
{}


// End of file
