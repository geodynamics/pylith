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

const double pylith::feassemble::ElasticityExplicitData2DQuadratic::_gravityVec[] = {
  0.00000000e+00, -1.00000000e+08,
};

const double pylith::feassemble::ElasticityExplicitData2DQuadratic::_vertices[] = {
 -1.00000000e+00, -1.00000000e+00,
  1.00000000e+00,  2.00000000e-01,
 -1.50000000e+00,  5.00000000e-01,
 -2.50000000e-01,  3.50000000e-01,
 -1.25000000e+00, -2.50000000e-01,
  0.00000000e+00, -4.00000000e-01,
};

const int pylith::feassemble::ElasticityExplicitData2DQuadratic::_cells[] = {
0,1,2,3,4,5,
};

const double pylith::feassemble::ElasticityExplicitData2DQuadratic::_verticesRef[] = {
 -1.00000000e+00, -1.00000000e+00,
  1.00000000e+00, -1.00000000e+00,
 -1.00000000e+00,  1.00000000e+00,
  0.00000000e+00,  0.00000000e+00,
 -1.00000000e+00,  0.00000000e+00,
  0.00000000e+00, -1.00000000e+00,
};

const double pylith::feassemble::ElasticityExplicitData2DQuadratic::_quadPts[] = {
  0.00000000e+00, -7.50000000e-01,
 -7.50000000e-01,  0.00000000e+00,
 -7.50000000e-01, -7.50000000e-01,
};

const double pylith::feassemble::ElasticityExplicitData2DQuadratic::_quadWts[] = {
  6.66666667e-01,  6.66666667e-01,  6.66666667e-01,
};

const double pylith::feassemble::ElasticityExplicitData2DQuadratic::_basis[] = {
 -9.37500000e-02,  0.00000000e+00,
 -9.37500000e-02,  2.50000000e-01,
  1.87500000e-01,  7.50000000e-01,
 -9.37500000e-02, -9.37500000e-02,
  0.00000000e+00,  2.50000000e-01,
  7.50000000e-01,  1.87500000e-01,
  3.75000000e-01, -9.37500000e-02,
 -9.37500000e-02,  6.25000000e-02,
  3.75000000e-01,  3.75000000e-01,
};

const double pylith::feassemble::ElasticityExplicitData2DQuadratic::_basisDerivRef[] = {
 -2.50000000e-01, -2.50000000e-01,
  5.00000000e-01,  0.00000000e+00,
  0.00000000e+00, -2.50000000e-01,
  2.50000000e-01,  1.00000000e+00,
 -2.50000000e-01,  5.00000000e-01,
 -2.50000000e-01, -1.00000000e+00,
 -2.50000000e-01, -2.50000000e-01,
 -2.50000000e-01,  0.00000000e+00,
  0.00000000e+00,  5.00000000e-01,
  1.00000000e+00,  2.50000000e-01,
 -1.00000000e+00, -2.50000000e-01,
  5.00000000e-01, -2.50000000e-01,
 -1.00000000e+00, -1.00000000e+00,
 -2.50000000e-01,  0.00000000e+00,
  0.00000000e+00, -2.50000000e-01,
  2.50000000e-01,  2.50000000e-01,
 -2.50000000e-01,  1.25000000e+00,
  1.25000000e+00, -2.50000000e-01,
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
  0.00000000e+00, -1.00000000e-01,
 -3.00000000e-01, -2.00000000e-01,
  1.00000000e-01,  4.00000000e-01,
 -2.00000000e-01, -6.00000000e-01,
};

const double pylith::feassemble::ElasticityExplicitData2DQuadratic::_valsResidual[] = {
 -5.90190847e+09,  4.71767000e+10,
 -1.39395410e+10, -6.36899121e+09,
  4.56827962e+08,  1.65078960e+10,
  2.42224578e+10, -8.30513314e+09,
 -2.84017435e+10, -6.85808242e+10,
  2.35658760e+10,  1.95703525e+10,
};

const double pylith::feassemble::ElasticityExplicitData2DQuadratic::_valsJacobian[] = {
  2.37304687e+06,  0.00000000e+00,
 -3.95507812e+05,  0.00000000e+00,
 -3.95507812e+05,  0.00000000e+00,
 -3.51562500e+05,  0.00000000e+00,
  7.91015625e+05,  0.00000000e+00,
  7.91015625e+05,  0.00000000e+00,
  0.00000000e+00,  2.37304687e+06,
  0.00000000e+00, -3.95507812e+05,
  0.00000000e+00, -3.95507812e+05,
  0.00000000e+00, -3.51562500e+05,
  0.00000000e+00,  7.91015625e+05,
  0.00000000e+00,  7.91015625e+05,
 -3.95507812e+05,  0.00000000e+00,
  2.63671875e+05,  0.00000000e+00,
  1.31835937e+05,  0.00000000e+00,
 -4.39453125e+05,  0.00000000e+00,
 -1.58203125e+06,  0.00000000e+00,
 -7.91015625e+05,  0.00000000e+00,
  0.00000000e+00, -3.95507812e+05,
  0.00000000e+00,  2.63671875e+05,
  0.00000000e+00,  1.31835937e+05,
  0.00000000e+00, -4.39453125e+05,
  0.00000000e+00, -1.58203125e+06,
  0.00000000e+00, -7.91015625e+05,
 -3.95507812e+05,  0.00000000e+00,
  1.31835937e+05,  0.00000000e+00,
  2.63671875e+05,  0.00000000e+00,
 -4.39453125e+05,  0.00000000e+00,
 -7.91015625e+05,  0.00000000e+00,
 -1.58203125e+06,  0.00000000e+00,
  0.00000000e+00, -3.95507812e+05,
  0.00000000e+00,  1.31835937e+05,
  0.00000000e+00,  2.63671875e+05,
  0.00000000e+00, -4.39453125e+05,
  0.00000000e+00, -7.91015625e+05,
  0.00000000e+00, -1.58203125e+06,
 -3.51562500e+05,  0.00000000e+00,
 -4.39453125e+05,  0.00000000e+00,
 -4.39453125e+05,  0.00000000e+00,
  1.93359375e+06,  0.00000000e+00,
  3.86718750e+06,  0.00000000e+00,
  3.86718750e+06,  0.00000000e+00,
  0.00000000e+00, -3.51562500e+05,
  0.00000000e+00, -4.39453125e+05,
  0.00000000e+00, -4.39453125e+05,
  0.00000000e+00,  1.93359375e+06,
  0.00000000e+00,  3.86718750e+06,
  0.00000000e+00,  3.86718750e+06,
  7.91015625e+05,  0.00000000e+00,
 -1.58203125e+06,  0.00000000e+00,
 -7.91015625e+05,  0.00000000e+00,
  3.86718750e+06,  0.00000000e+00,
  1.10742188e+07,  0.00000000e+00,
  6.32812500e+06,  0.00000000e+00,
  0.00000000e+00,  7.91015625e+05,
  0.00000000e+00, -1.58203125e+06,
  0.00000000e+00, -7.91015625e+05,
  0.00000000e+00,  3.86718750e+06,
  0.00000000e+00,  1.10742188e+07,
  0.00000000e+00,  6.32812500e+06,
  7.91015625e+05,  0.00000000e+00,
 -7.91015625e+05,  0.00000000e+00,
 -1.58203125e+06,  0.00000000e+00,
  3.86718750e+06,  0.00000000e+00,
  6.32812500e+06,  0.00000000e+00,
  1.10742188e+07,  0.00000000e+00,
  0.00000000e+00,  7.91015625e+05,
  0.00000000e+00, -7.91015625e+05,
  0.00000000e+00, -1.58203125e+06,
  0.00000000e+00,  3.86718750e+06,
  0.00000000e+00,  6.32812500e+06,
  0.00000000e+00,  1.10742188e+07,
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
  gravityVec = const_cast<double*>(_gravityVec);
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

pylith::feassemble::ElasticityExplicitData2DQuadratic::~ElasticityExplicitData2DQuadratic(void)
{}


// End of file
