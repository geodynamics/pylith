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

const double pylith::feassemble::ElasticityImplicitData2DQuadratic::_verticesRef[] = {
 -1.00000000e+00, -1.00000000e+00,
  1.00000000e+00, -1.00000000e+00,
 -1.00000000e+00,  1.00000000e+00,
  0.00000000e+00,  0.00000000e+00,
 -1.00000000e+00,  0.00000000e+00,
  0.00000000e+00, -1.00000000e+00,
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

const double pylith::feassemble::ElasticityImplicitData2DQuadratic::_basisDerivRef[] = {
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
 -7.03759529e+08, -4.10320114e+10,
 -3.70610450e+10,  2.11260736e+10,
 -1.31513835e+10, -1.25323831e+11,
  4.33964587e+10,  1.47829671e+11,
 -5.10800854e+10, -2.16829139e+11,
  5.85998147e+10,  2.14229238e+11,
};

const double pylith::feassemble::ElasticityImplicitData2DQuadratic::_valsJacobian[] = {
  3.10620755e+10,  1.51816443e+10,
  4.03724615e+09,  4.59327234e+09,
  7.44201344e+09, -1.16900072e+09,
 -1.88162512e+10, -1.81941422e+10,
 -5.31233692e+09, -2.78656421e+09,
 -1.84127470e+10,  2.37479043e+09,
  1.51816443e+10,  7.32318930e+10,
  4.59327234e+09, -7.39325357e+08,
 -1.16900072e+09,  2.72257870e+10,
 -1.81941422e+10, -4.08528660e+09,
 -2.78656421e+09, -6.82095745e+09,
  2.37479043e+09, -8.88121106e+10,
  4.03724615e+09,  4.59327234e+09,
  5.52513202e+10, -1.36681382e+10,
  1.16542644e+10, -7.01226095e+09,
  3.42514288e+09, -2.14282837e+10,
 -6.14617630e+10,  3.14248548e+10,
 -1.29062107e+10,  6.09055570e+09,
  4.59327234e+09, -7.39325357e+08,
 -1.36681382e+10,  2.59749439e+10,
 -7.01226095e+09,  7.52420433e+09,
 -2.14282837e+10,  2.63944418e+10,
  3.14248548e+10, -5.19138623e+10,
  6.09055570e+09, -7.24040229e+09,
  7.44201344e+09, -1.16900072e+09,
  1.16542644e+10, -7.01226095e+09,
  4.11450311e+10, -2.14269305e+10,
 -4.46688705e+09, -7.49680571e+07,
 -2.01418211e+10,  2.05378765e+10,
 -3.56326007e+10,  9.14528366e+09,
 -1.16900072e+09,  2.72257870e+10,
 -7.01226095e+09,  7.52420433e+09,
 -2.14269305e+10,  9.10320049e+10,
 -7.49680571e+07, -1.35889491e+10,
  2.05378765e+10,  6.67218742e+09,
  9.14528366e+09, -1.18865235e+11,
 -1.88162512e+10, -1.81941422e+10,
  3.42514288e+09, -2.14282837e+10,
 -4.46688705e+09, -7.49680571e+07,
  9.34772297e+10, -1.75952640e+10,
 -3.26438450e+10,  2.31583412e+10,
 -4.09753893e+10,  3.41343167e+10,
 -1.81941422e+10, -4.08528660e+09,
 -2.14282837e+10,  2.63944418e+10,
 -7.49680571e+07, -1.35889491e+10,
 -1.75952640e+10,  1.28993398e+11,
  2.31583412e+10, -1.32670104e+11,
  3.41343167e+10, -5.04349923e+09,
 -5.31233692e+09, -2.78656421e+09,
 -6.14617630e+10,  3.14248548e+10,
 -2.01418211e+10,  2.05378765e+10,
 -3.26438450e+10,  2.31583412e+10,
  1.28924268e+11, -2.92151596e+10,
 -9.36450187e+09, -4.31193486e+10,
 -2.78656421e+09, -6.82095745e+09,
  3.14248548e+10, -5.19138623e+10,
  2.05378765e+10,  6.67218742e+09,
  2.31583412e+10, -1.32670104e+11,
 -2.92151596e+10,  1.94678071e+11,
 -4.31193486e+10, -9.94533453e+09,
 -1.84127470e+10,  2.37479043e+09,
 -1.29062107e+10,  6.09055570e+09,
 -3.56326007e+10,  9.14528366e+09,
 -4.09753893e+10,  3.41343167e+10,
 -9.36450187e+09, -4.31193486e+10,
  1.17291450e+11, -8.62559790e+09,
  2.37479043e+09, -8.88121106e+10,
  6.09055570e+09, -7.24040229e+09,
  9.14528366e+09, -1.18865235e+11,
  3.41343167e+10, -5.04349923e+09,
 -4.31193486e+10, -9.94533453e+09,
 -8.62559790e+09,  2.29906581e+11,
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

pylith::feassemble::ElasticityImplicitData2DQuadratic::~ElasticityImplicitData2DQuadratic(void)
{}


// End of file
