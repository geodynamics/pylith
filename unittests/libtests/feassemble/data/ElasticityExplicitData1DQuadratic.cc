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

#include "ElasticityExplicitData1DQuadratic.hh"

const int pylith::feassemble::ElasticityExplicitData1DQuadratic::_spaceDim = 1;

const int pylith::feassemble::ElasticityExplicitData1DQuadratic::_cellDim = 1;

const int pylith::feassemble::ElasticityExplicitData1DQuadratic::_numVertices = 3;

const int pylith::feassemble::ElasticityExplicitData1DQuadratic::_numCells = 1;

const int pylith::feassemble::ElasticityExplicitData1DQuadratic::_numBasis = 3;

const int pylith::feassemble::ElasticityExplicitData1DQuadratic::_numQuadPts = 2;

const char* pylith::feassemble::ElasticityExplicitData1DQuadratic::_matType = "ElasticStrain1D";

const char* pylith::feassemble::ElasticityExplicitData1DQuadratic::_matDBFilename = "data/elasticstrain1d.spatialdb";

const int pylith::feassemble::ElasticityExplicitData1DQuadratic::_matId = 0;

const char* pylith::feassemble::ElasticityExplicitData1DQuadratic::_matLabel = "elastic strain 1-D";

const double pylith::feassemble::ElasticityExplicitData1DQuadratic::_dt =   1.00000000e-02;

const double pylith::feassemble::ElasticityExplicitData1DQuadratic::_vertices[] = {
 -2.50000000e-01,
  8.75000000e-01,
  2.00000000e+00,
};

const int pylith::feassemble::ElasticityExplicitData1DQuadratic::_cells[] = {
0,1,2,
};

const double pylith::feassemble::ElasticityExplicitData1DQuadratic::_verticesRef[] = {
 -1.00000000e+00,
  1.00000000e+00,
  0.00000000e+00,
};

const double pylith::feassemble::ElasticityExplicitData1DQuadratic::_quadPts[] = {
 -5.77350269e-01,
  5.77350269e-01,
};

const double pylith::feassemble::ElasticityExplicitData1DQuadratic::_quadWts[] = {
  1.00000000e+00,  1.00000000e+00,
};

const double pylith::feassemble::ElasticityExplicitData1DQuadratic::_basis[] = {
  4.55341801e-01,
  6.66666667e-01,
 -1.22008468e-01,
 -1.22008468e-01,
  6.66666667e-01,
  4.55341801e-01,
};

const double pylith::feassemble::ElasticityExplicitData1DQuadratic::_basisDerivRef[] = {
 -1.07735027e+00,
  1.15470054e+00,
 -7.73502692e-02,
  7.73502692e-02,
 -1.15470054e+00,
  1.07735027e+00,
};

const double pylith::feassemble::ElasticityExplicitData1DQuadratic::_fieldTpdt[] = {
  1.20000000e+00,
  0.00000000e+00,
  1.70000000e+00,
};

const double pylith::feassemble::ElasticityExplicitData1DQuadratic::_fieldT[] = {
  1.10000000e+00,
  1.00000000e-01,
  1.50000000e+00,
};

const double pylith::feassemble::ElasticityExplicitData1DQuadratic::_fieldTmdt[] = {
  1.00000000e+00,
  1.00000000e-01,
  1.30000000e+00,
};

const double pylith::feassemble::ElasticityExplicitData1DQuadratic::_valsResidual[] = {
 -1.11997188e+11,
  2.56020625e+11,
 -1.43992500e+11,
};

const double pylith::feassemble::ElasticityExplicitData1DQuadratic::_valsJacobian[] = {
  6.25000000e+06,
  6.25000000e+06,
 -3.12500000e+06,
  6.25000000e+06,
  2.50000000e+07,
  6.25000000e+06,
 -3.12500000e+06,
  6.25000000e+06,
  6.25000000e+06,
};

pylith::feassemble::ElasticityExplicitData1DQuadratic::ElasticityExplicitData1DQuadratic(void)
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

pylith::feassemble::ElasticityExplicitData1DQuadratic::~ElasticityExplicitData1DQuadratic(void)
{}


// End of file
