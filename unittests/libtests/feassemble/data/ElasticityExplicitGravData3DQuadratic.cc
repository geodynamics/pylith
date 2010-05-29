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
// This file was generated from python application elasticityexplicitapp.

#include "ElasticityExplicitGravData3DQuadratic.hh"

const int pylith::feassemble::ElasticityExplicitGravData3DQuadratic::_spaceDim = 3;

const int pylith::feassemble::ElasticityExplicitGravData3DQuadratic::_cellDim = 3;

const int pylith::feassemble::ElasticityExplicitGravData3DQuadratic::_numVertices = 10;

const int pylith::feassemble::ElasticityExplicitGravData3DQuadratic::_numCells = 1;

const int pylith::feassemble::ElasticityExplicitGravData3DQuadratic::_numBasis = 10;

const int pylith::feassemble::ElasticityExplicitGravData3DQuadratic::_numQuadPts = 4;

const char* pylith::feassemble::ElasticityExplicitGravData3DQuadratic::_matType = "ElasticIsotropic3D";

const char* pylith::feassemble::ElasticityExplicitGravData3DQuadratic::_matDBFilename = "data/elasticisotropic3d.spatialdb";

const int pylith::feassemble::ElasticityExplicitGravData3DQuadratic::_matId = 0;

const char* pylith::feassemble::ElasticityExplicitGravData3DQuadratic::_matLabel = "elastic isotropic 3-D";

const double pylith::feassemble::ElasticityExplicitGravData3DQuadratic::_dt =   1.00000000e-02;

const double pylith::feassemble::ElasticityExplicitGravData3DQuadratic::_gravityVec[] = {
  0.00000000e+00,  0.00000000e+00, -1.00000000e+08,
};

const double pylith::feassemble::ElasticityExplicitGravData3DQuadratic::_vertices[] = {
 -5.00000000e-01, -2.00000000e+00, -1.00000000e+00,
  2.00000000e+00, -2.00000000e+00, -5.00000000e-01,
  1.00000000e+00,  1.00000000e+00,  0.00000000e+00,
  2.00000000e-01,  5.00000000e-01,  2.00000000e+00,
  1.50000000e+00, -5.00000000e-01, -2.50000000e-01,
  2.50000000e-01, -5.00000000e-01, -5.00000000e-01,
  7.50000000e-01, -2.00000000e+00, -7.50000000e-01,
 -1.50000000e-01, -7.50000000e-01,  5.00000000e-01,
  1.10000000e+00, -7.50000000e-01,  7.50000000e-01,
  6.00000000e-01,  7.50000000e-01,  1.00000000e+00,
};

const int pylith::feassemble::ElasticityExplicitGravData3DQuadratic::_cells[] = {
0,1,2,3,4,5,6,7,8,9,
};

const double pylith::feassemble::ElasticityExplicitGravData3DQuadratic::_verticesRef[] = {
 -1.00000000e+00, -1.00000000e+00, -1.00000000e+00,
  1.00000000e+00, -1.00000000e+00, -1.00000000e+00,
 -1.00000000e+00,  1.00000000e+00, -1.00000000e+00,
 -1.00000000e+00, -1.00000000e+00,  1.00000000e+00,
  0.00000000e+00,  0.00000000e+00, -1.00000000e+00,
 -1.00000000e+00,  0.00000000e+00, -1.00000000e+00,
  0.00000000e+00, -1.00000000e+00, -1.00000000e+00,
 -1.00000000e+00, -1.00000000e+00,  0.00000000e+00,
  0.00000000e+00, -1.00000000e+00,  0.00000000e+00,
 -1.00000000e+00,  0.00000000e+00,  0.00000000e+00,
};

const double pylith::feassemble::ElasticityExplicitGravData3DQuadratic::_quadPts[] = {
 -8.00000000e-01, -8.00000000e-01, -8.00000000e-01,
  5.00000000e-01, -8.00000000e-01, -8.00000000e-01,
 -8.00000000e-01,  5.00000000e-01, -8.00000000e-01,
 -8.00000000e-01, -8.00000000e-01,  5.00000000e-01,
};

const double pylith::feassemble::ElasticityExplicitGravData3DQuadratic::_quadWts[] = {
  3.33333333e-01,  3.33333333e-01,  3.33333333e-01,  3.33333333e-01,
};

const double pylith::feassemble::ElasticityExplicitGravData3DQuadratic::_basis[] = {
  2.80000000e-01, -8.00000000e-02, -8.00000000e-02,
 -8.00000000e-02,  4.00000000e-02,  2.80000000e-01,
  2.80000000e-01,  2.80000000e-01,  4.00000000e-02,
  4.00000000e-02, -4.50000000e-02,  3.75000000e-01,
 -8.00000000e-02, -8.00000000e-02,  3.00000000e-01,
  2.00000000e-02,  1.50000000e-01,  2.00000000e-02,
  3.00000000e-01,  4.00000000e-02, -4.50000000e-02,
 -8.00000000e-02,  3.75000000e-01, -8.00000000e-02,
  3.00000000e-01,  1.50000000e-01,  2.00000000e-02,
  2.00000000e-02,  4.00000000e-02,  3.00000000e-01,
 -4.50000000e-02, -8.00000000e-02, -8.00000000e-02,
  3.75000000e-01,  4.00000000e-02,  2.00000000e-02,
  2.00000000e-02,  1.50000000e-01,  3.00000000e-01,
  3.00000000e-01,};

const double pylith::feassemble::ElasticityExplicitGravData3DQuadratic::_basisDerivRef[] = {
 -9.00000000e-01, -9.00000000e-01, -9.00000000e-01,
 -3.00000000e-01,  0.00000000e+00,  0.00000000e+00,
  0.00000000e+00, -3.00000000e-01,  0.00000000e+00,
  0.00000000e+00,  0.00000000e+00, -3.00000000e-01,
  2.00000000e-01,  2.00000000e-01,  0.00000000e+00,
 -2.00000000e-01,  1.20000000e+00, -2.00000000e-01,
  1.20000000e+00, -2.00000000e-01, -2.00000000e-01,
 -2.00000000e-01, -2.00000000e-01,  1.20000000e+00,
  2.00000000e-01,  0.00000000e+00,  2.00000000e-01,
  0.00000000e+00,  2.00000000e-01,  2.00000000e-01,
  4.00000000e-01,  4.00000000e-01,  4.00000000e-01,
  1.00000000e+00,  0.00000000e+00,  0.00000000e+00,
  0.00000000e+00, -3.00000000e-01,  0.00000000e+00,
  0.00000000e+00,  0.00000000e+00, -3.00000000e-01,
  2.00000000e-01,  1.50000000e+00,  0.00000000e+00,
 -2.00000000e-01, -1.00000000e-01, -2.00000000e-01,
 -1.40000000e+00, -1.50000000e+00, -1.50000000e+00,
 -2.00000000e-01, -2.00000000e-01, -1.00000000e-01,
  2.00000000e-01,  0.00000000e+00,  1.50000000e+00,
  0.00000000e+00,  2.00000000e-01,  2.00000000e-01,
  4.00000000e-01,  4.00000000e-01,  4.00000000e-01,
 -3.00000000e-01,  0.00000000e+00,  0.00000000e+00,
  0.00000000e+00,  1.00000000e+00,  0.00000000e+00,
  0.00000000e+00,  0.00000000e+00, -3.00000000e-01,
  1.50000000e+00,  2.00000000e-01,  0.00000000e+00,
 -1.50000000e+00, -1.40000000e+00, -1.50000000e+00,
 -1.00000000e-01, -2.00000000e-01, -2.00000000e-01,
 -2.00000000e-01, -2.00000000e-01, -1.00000000e-01,
  2.00000000e-01,  0.00000000e+00,  2.00000000e-01,
  0.00000000e+00,  2.00000000e-01,  1.50000000e+00,
  4.00000000e-01,  4.00000000e-01,  4.00000000e-01,
 -3.00000000e-01,  0.00000000e+00,  0.00000000e+00,
  0.00000000e+00, -3.00000000e-01,  0.00000000e+00,
  0.00000000e+00,  0.00000000e+00,  1.00000000e+00,
  2.00000000e-01,  2.00000000e-01,  0.00000000e+00,
 -2.00000000e-01, -1.00000000e-01, -2.00000000e-01,
 -1.00000000e-01, -2.00000000e-01, -2.00000000e-01,
 -1.50000000e+00, -1.50000000e+00, -1.40000000e+00,
  1.50000000e+00,  0.00000000e+00,  2.00000000e-01,
  0.00000000e+00,  1.50000000e+00,  2.00000000e-01,
};

const double pylith::feassemble::ElasticityExplicitGravData3DQuadratic::_fieldTIncr[] = {
  3.00000000e-01, -4.00000000e-01, -4.00000000e-01,
 -6.00000000e-01,  8.00000000e-01,  2.00000000e-01,
  5.00000000e-01,  5.00000000e-01,  7.00000000e-01,
 -7.00000000e-01, -5.00000000e-01, -7.00000000e-01,
 -6.00000000e-01, -3.00000000e-01,  8.00000000e-01,
 -4.00000000e-01, -8.00000000e-01, -5.00000000e-01,
  7.00000000e-01,  8.00000000e-01, -5.00000000e-01,
 -5.00000000e-01, -5.00000000e-01, -7.00000000e-01,
 -3.00000000e-01, -9.00000000e-01,  8.00000000e-01,
 -1.00000000e-01,  5.00000000e-01, -9.00000000e-01,
};

const double pylith::feassemble::ElasticityExplicitGravData3DQuadratic::_fieldT[] = {
  1.00000000e-01, -2.00000000e-01, -6.00000000e-01,
 -3.00000000e-01,  4.00000000e-01,  9.00000000e-01,
  6.00000000e-01,  8.00000000e-01,  5.00000000e-01,
 -8.00000000e-01, -6.00000000e-01, -8.00000000e-01,
 -0.00000000e+00, -2.00000000e-01,  6.00000000e-01,
 -4.00000000e-01, -7.00000000e-01, -2.00000000e-01,
  7.00000000e-01,  6.00000000e-01, -1.00000000e-01,
 -4.00000000e-01, -3.00000000e-01, -3.00000000e-01,
 -7.00000000e-01, -6.00000000e-01,  1.00000000e-01,
 -9.00000000e-01,  3.00000000e-01, -8.00000000e-01,
};

const double pylith::feassemble::ElasticityExplicitGravData3DQuadratic::_fieldTmdt[] = {
  2.00000000e-01, -3.00000000e-01, -1.00000000e-01,
 -4.00000000e-01,  2.00000000e-01,  3.00000000e-01,
 -5.00000000e-01,  2.00000000e-01,  2.00000000e-01,
 -3.00000000e-01, -8.00000000e-01, -3.00000000e-01,
 -5.00000000e-01, -9.00000000e-01,  4.00000000e-01,
 -3.00000000e-01, -6.00000000e-01, -8.00000000e-01,
  9.00000000e-01,  5.00000000e-01, -2.00000000e-01,
 -7.00000000e-01, -2.00000000e-01, -9.00000000e-01,
 -5.00000000e-01, -8.00000000e-01,  4.00000000e-01,
 -4.00000000e-01,  5.00000000e-01, -7.00000000e-01,
};

const double pylith::feassemble::ElasticityExplicitGravData3DQuadratic::_valsResidual[] = {
  2.17437265e+10, -9.08080342e+09, -5.09278842e+09,
  6.22359327e+10, -2.47895200e+10, -1.65925044e+10,
 -4.60088086e+10, -5.34117841e+10, -5.55302367e+10,
 -1.00528639e+10,  4.19879745e+10, -6.27435491e+10,
 -8.42948094e+09,  6.33090926e+10, -1.24369193e+11,
  6.62259132e+10,  1.09240806e+11, -1.20874795e+11,
 -9.59434641e+10, -5.29480086e+10, -1.08714691e+11,
 -1.25307052e+10, -6.06736178e+10, -7.32224427e+10,
 -3.24399696e+10,  5.14120021e+10, -1.13619372e+11,
  5.52136379e+10, -6.50240950e+10, -3.06843461e+10,
};

const double pylith::feassemble::ElasticityExplicitGravData3DQuadratic::_valsJacobian[] = {
  1.50251107e+06,  0.00000000e+00,  0.00000000e+00,
 -5.70500651e+05,  0.00000000e+00,  0.00000000e+00,
 -5.70500651e+05,  0.00000000e+00,  0.00000000e+00,
 -5.70500651e+05,  0.00000000e+00,  0.00000000e+00,
 -3.13041667e+05,  0.00000000e+00,  0.00000000e+00,
  1.24238411e+06,  0.00000000e+00,  0.00000000e+00,
  1.24238411e+06,  0.00000000e+00,  0.00000000e+00,
  1.24238411e+06,  0.00000000e+00,  0.00000000e+00,
 -3.13041667e+05,  0.00000000e+00,  0.00000000e+00,
 -3.13041667e+05,  0.00000000e+00,  0.00000000e+00,
  0.00000000e+00,  1.50251107e+06,  0.00000000e+00,
  0.00000000e+00, -5.70500651e+05,  0.00000000e+00,
  0.00000000e+00, -5.70500651e+05,  0.00000000e+00,
  0.00000000e+00, -5.70500651e+05,  0.00000000e+00,
  0.00000000e+00, -3.13041667e+05,  0.00000000e+00,
  0.00000000e+00,  1.24238411e+06,  0.00000000e+00,
  0.00000000e+00,  1.24238411e+06,  0.00000000e+00,
  0.00000000e+00,  1.24238411e+06,  0.00000000e+00,
  0.00000000e+00, -3.13041667e+05,  0.00000000e+00,
  0.00000000e+00, -3.13041667e+05,  0.00000000e+00,
  0.00000000e+00,  0.00000000e+00,  1.50251107e+06,
  0.00000000e+00,  0.00000000e+00, -5.70500651e+05,
  0.00000000e+00,  0.00000000e+00, -5.70500651e+05,
  0.00000000e+00,  0.00000000e+00, -5.70500651e+05,
  0.00000000e+00,  0.00000000e+00, -3.13041667e+05,
  0.00000000e+00,  0.00000000e+00,  1.24238411e+06,
  0.00000000e+00,  0.00000000e+00,  1.24238411e+06,
  0.00000000e+00,  0.00000000e+00,  1.24238411e+06,
  0.00000000e+00,  0.00000000e+00, -3.13041667e+05,
  0.00000000e+00,  0.00000000e+00, -3.13041667e+05,
 -5.70500651e+05,  0.00000000e+00,  0.00000000e+00,
  2.84272070e+06,  0.00000000e+00,  0.00000000e+00,
 -8.39520833e+05,  0.00000000e+00,  0.00000000e+00,
 -8.39520833e+05,  0.00000000e+00,  0.00000000e+00,
  1.46026823e+06,  0.00000000e+00,  0.00000000e+00,
 -5.06914062e+05,  0.00000000e+00,  0.00000000e+00,
  5.45154948e+05,  0.00000000e+00,  0.00000000e+00,
 -5.06914062e+05,  0.00000000e+00,  0.00000000e+00,
  1.46026823e+06,  0.00000000e+00,  0.00000000e+00,
 -6.43869792e+05,  0.00000000e+00,  0.00000000e+00,
  0.00000000e+00, -5.70500651e+05,  0.00000000e+00,
  0.00000000e+00,  2.84272070e+06,  0.00000000e+00,
  0.00000000e+00, -8.39520833e+05,  0.00000000e+00,
  0.00000000e+00, -8.39520833e+05,  0.00000000e+00,
  0.00000000e+00,  1.46026823e+06,  0.00000000e+00,
  0.00000000e+00, -5.06914062e+05,  0.00000000e+00,
  0.00000000e+00,  5.45154948e+05,  0.00000000e+00,
  0.00000000e+00, -5.06914062e+05,  0.00000000e+00,
  0.00000000e+00,  1.46026823e+06,  0.00000000e+00,
  0.00000000e+00, -6.43869792e+05,  0.00000000e+00,
  0.00000000e+00,  0.00000000e+00, -5.70500651e+05,
  0.00000000e+00,  0.00000000e+00,  2.84272070e+06,
  0.00000000e+00,  0.00000000e+00, -8.39520833e+05,
  0.00000000e+00,  0.00000000e+00, -8.39520833e+05,
  0.00000000e+00,  0.00000000e+00,  1.46026823e+06,
  0.00000000e+00,  0.00000000e+00, -5.06914062e+05,
  0.00000000e+00,  0.00000000e+00,  5.45154948e+05,
  0.00000000e+00,  0.00000000e+00, -5.06914062e+05,
  0.00000000e+00,  0.00000000e+00,  1.46026823e+06,
  0.00000000e+00,  0.00000000e+00, -6.43869792e+05,
 -5.70500651e+05,  0.00000000e+00,  0.00000000e+00,
 -8.39520833e+05,  0.00000000e+00,  0.00000000e+00,
  2.84272070e+06,  0.00000000e+00,  0.00000000e+00,
 -8.39520833e+05,  0.00000000e+00,  0.00000000e+00,
  1.46026823e+06,  0.00000000e+00,  0.00000000e+00,
  5.45154948e+05,  0.00000000e+00,  0.00000000e+00,
 -5.06914062e+05,  0.00000000e+00,  0.00000000e+00,
 -5.06914062e+05,  0.00000000e+00,  0.00000000e+00,
 -6.43869792e+05,  0.00000000e+00,  0.00000000e+00,
  1.46026823e+06,  0.00000000e+00,  0.00000000e+00,
  0.00000000e+00, -5.70500651e+05,  0.00000000e+00,
  0.00000000e+00, -8.39520833e+05,  0.00000000e+00,
  0.00000000e+00,  2.84272070e+06,  0.00000000e+00,
  0.00000000e+00, -8.39520833e+05,  0.00000000e+00,
  0.00000000e+00,  1.46026823e+06,  0.00000000e+00,
  0.00000000e+00,  5.45154948e+05,  0.00000000e+00,
  0.00000000e+00, -5.06914062e+05,  0.00000000e+00,
  0.00000000e+00, -5.06914062e+05,  0.00000000e+00,
  0.00000000e+00, -6.43869792e+05,  0.00000000e+00,
  0.00000000e+00,  1.46026823e+06,  0.00000000e+00,
  0.00000000e+00,  0.00000000e+00, -5.70500651e+05,
  0.00000000e+00,  0.00000000e+00, -8.39520833e+05,
  0.00000000e+00,  0.00000000e+00,  2.84272070e+06,
  0.00000000e+00,  0.00000000e+00, -8.39520833e+05,
  0.00000000e+00,  0.00000000e+00,  1.46026823e+06,
  0.00000000e+00,  0.00000000e+00,  5.45154948e+05,
  0.00000000e+00,  0.00000000e+00, -5.06914062e+05,
  0.00000000e+00,  0.00000000e+00, -5.06914062e+05,
  0.00000000e+00,  0.00000000e+00, -6.43869792e+05,
  0.00000000e+00,  0.00000000e+00,  1.46026823e+06,
 -5.70500651e+05,  0.00000000e+00,  0.00000000e+00,
 -8.39520833e+05,  0.00000000e+00,  0.00000000e+00,
 -8.39520833e+05,  0.00000000e+00,  0.00000000e+00,
  2.84272070e+06,  0.00000000e+00,  0.00000000e+00,
 -6.43869792e+05,  0.00000000e+00,  0.00000000e+00,
 -5.06914062e+05,  0.00000000e+00,  0.00000000e+00,
 -5.06914062e+05,  0.00000000e+00,  0.00000000e+00,
  5.45154948e+05,  0.00000000e+00,  0.00000000e+00,
  1.46026823e+06,  0.00000000e+00,  0.00000000e+00,
  1.46026823e+06,  0.00000000e+00,  0.00000000e+00,
  0.00000000e+00, -5.70500651e+05,  0.00000000e+00,
  0.00000000e+00, -8.39520833e+05,  0.00000000e+00,
  0.00000000e+00, -8.39520833e+05,  0.00000000e+00,
  0.00000000e+00,  2.84272070e+06,  0.00000000e+00,
  0.00000000e+00, -6.43869792e+05,  0.00000000e+00,
  0.00000000e+00, -5.06914062e+05,  0.00000000e+00,
  0.00000000e+00, -5.06914062e+05,  0.00000000e+00,
  0.00000000e+00,  5.45154948e+05,  0.00000000e+00,
  0.00000000e+00,  1.46026823e+06,  0.00000000e+00,
  0.00000000e+00,  1.46026823e+06,  0.00000000e+00,
  0.00000000e+00,  0.00000000e+00, -5.70500651e+05,
  0.00000000e+00,  0.00000000e+00, -8.39520833e+05,
  0.00000000e+00,  0.00000000e+00, -8.39520833e+05,
  0.00000000e+00,  0.00000000e+00,  2.84272070e+06,
  0.00000000e+00,  0.00000000e+00, -6.43869792e+05,
  0.00000000e+00,  0.00000000e+00, -5.06914062e+05,
  0.00000000e+00,  0.00000000e+00, -5.06914062e+05,
  0.00000000e+00,  0.00000000e+00,  5.45154948e+05,
  0.00000000e+00,  0.00000000e+00,  1.46026823e+06,
  0.00000000e+00,  0.00000000e+00,  1.46026823e+06,
 -3.13041667e+05,  0.00000000e+00,  0.00000000e+00,
  1.46026823e+06,  0.00000000e+00,  0.00000000e+00,
  1.46026823e+06,  0.00000000e+00,  0.00000000e+00,
 -6.43869792e+05,  0.00000000e+00,  0.00000000e+00,
  3.25847917e+06,  0.00000000e+00,  0.00000000e+00,
  1.12054688e+06,  0.00000000e+00,  0.00000000e+00,
  1.12054688e+06,  0.00000000e+00,  0.00000000e+00,
  5.19364583e+05,  0.00000000e+00,  0.00000000e+00,
  2.05611458e+06,  0.00000000e+00,  0.00000000e+00,
  2.05611458e+06,  0.00000000e+00,  0.00000000e+00,
  0.00000000e+00, -3.13041667e+05,  0.00000000e+00,
  0.00000000e+00,  1.46026823e+06,  0.00000000e+00,
  0.00000000e+00,  1.46026823e+06,  0.00000000e+00,
  0.00000000e+00, -6.43869792e+05,  0.00000000e+00,
  0.00000000e+00,  3.25847917e+06,  0.00000000e+00,
  0.00000000e+00,  1.12054688e+06,  0.00000000e+00,
  0.00000000e+00,  1.12054688e+06,  0.00000000e+00,
  0.00000000e+00,  5.19364583e+05,  0.00000000e+00,
  0.00000000e+00,  2.05611458e+06,  0.00000000e+00,
  0.00000000e+00,  2.05611458e+06,  0.00000000e+00,
  0.00000000e+00,  0.00000000e+00, -3.13041667e+05,
  0.00000000e+00,  0.00000000e+00,  1.46026823e+06,
  0.00000000e+00,  0.00000000e+00,  1.46026823e+06,
  0.00000000e+00,  0.00000000e+00, -6.43869792e+05,
  0.00000000e+00,  0.00000000e+00,  3.25847917e+06,
  0.00000000e+00,  0.00000000e+00,  1.12054688e+06,
  0.00000000e+00,  0.00000000e+00,  1.12054688e+06,
  0.00000000e+00,  0.00000000e+00,  5.19364583e+05,
  0.00000000e+00,  0.00000000e+00,  2.05611458e+06,
  0.00000000e+00,  0.00000000e+00,  2.05611458e+06,
  1.24238411e+06,  0.00000000e+00,  0.00000000e+00,
 -5.06914062e+05,  0.00000000e+00,  0.00000000e+00,
  5.45154948e+05,  0.00000000e+00,  0.00000000e+00,
 -5.06914062e+05,  0.00000000e+00,  0.00000000e+00,
  1.12054688e+06,  0.00000000e+00,  0.00000000e+00,
  1.80888281e+06,  0.00000000e+00,  0.00000000e+00,
  1.50829167e+06,  0.00000000e+00,  0.00000000e+00,
  1.50829167e+06,  0.00000000e+00,  0.00000000e+00,
  5.19364583e+05,  0.00000000e+00,  0.00000000e+00,
  1.12054688e+06,  0.00000000e+00,  0.00000000e+00,
  0.00000000e+00,  1.24238411e+06,  0.00000000e+00,
  0.00000000e+00, -5.06914062e+05,  0.00000000e+00,
  0.00000000e+00,  5.45154948e+05,  0.00000000e+00,
  0.00000000e+00, -5.06914062e+05,  0.00000000e+00,
  0.00000000e+00,  1.12054688e+06,  0.00000000e+00,
  0.00000000e+00,  1.80888281e+06,  0.00000000e+00,
  0.00000000e+00,  1.50829167e+06,  0.00000000e+00,
  0.00000000e+00,  1.50829167e+06,  0.00000000e+00,
  0.00000000e+00,  5.19364583e+05,  0.00000000e+00,
  0.00000000e+00,  1.12054688e+06,  0.00000000e+00,
  0.00000000e+00,  0.00000000e+00,  1.24238411e+06,
  0.00000000e+00,  0.00000000e+00, -5.06914062e+05,
  0.00000000e+00,  0.00000000e+00,  5.45154948e+05,
  0.00000000e+00,  0.00000000e+00, -5.06914062e+05,
  0.00000000e+00,  0.00000000e+00,  1.12054688e+06,
  0.00000000e+00,  0.00000000e+00,  1.80888281e+06,
  0.00000000e+00,  0.00000000e+00,  1.50829167e+06,
  0.00000000e+00,  0.00000000e+00,  1.50829167e+06,
  0.00000000e+00,  0.00000000e+00,  5.19364583e+05,
  0.00000000e+00,  0.00000000e+00,  1.12054688e+06,
  1.24238411e+06,  0.00000000e+00,  0.00000000e+00,
  5.45154948e+05,  0.00000000e+00,  0.00000000e+00,
 -5.06914062e+05,  0.00000000e+00,  0.00000000e+00,
 -5.06914062e+05,  0.00000000e+00,  0.00000000e+00,
  1.12054688e+06,  0.00000000e+00,  0.00000000e+00,
  1.50829167e+06,  0.00000000e+00,  0.00000000e+00,
  1.80888281e+06,  0.00000000e+00,  0.00000000e+00,
  1.50829167e+06,  0.00000000e+00,  0.00000000e+00,
  1.12054688e+06,  0.00000000e+00,  0.00000000e+00,
  5.19364583e+05,  0.00000000e+00,  0.00000000e+00,
  0.00000000e+00,  1.24238411e+06,  0.00000000e+00,
  0.00000000e+00,  5.45154948e+05,  0.00000000e+00,
  0.00000000e+00, -5.06914062e+05,  0.00000000e+00,
  0.00000000e+00, -5.06914062e+05,  0.00000000e+00,
  0.00000000e+00,  1.12054688e+06,  0.00000000e+00,
  0.00000000e+00,  1.50829167e+06,  0.00000000e+00,
  0.00000000e+00,  1.80888281e+06,  0.00000000e+00,
  0.00000000e+00,  1.50829167e+06,  0.00000000e+00,
  0.00000000e+00,  1.12054688e+06,  0.00000000e+00,
  0.00000000e+00,  5.19364583e+05,  0.00000000e+00,
  0.00000000e+00,  0.00000000e+00,  1.24238411e+06,
  0.00000000e+00,  0.00000000e+00,  5.45154948e+05,
  0.00000000e+00,  0.00000000e+00, -5.06914062e+05,
  0.00000000e+00,  0.00000000e+00, -5.06914062e+05,
  0.00000000e+00,  0.00000000e+00,  1.12054688e+06,
  0.00000000e+00,  0.00000000e+00,  1.50829167e+06,
  0.00000000e+00,  0.00000000e+00,  1.80888281e+06,
  0.00000000e+00,  0.00000000e+00,  1.50829167e+06,
  0.00000000e+00,  0.00000000e+00,  1.12054688e+06,
  0.00000000e+00,  0.00000000e+00,  5.19364583e+05,
  1.24238411e+06,  0.00000000e+00,  0.00000000e+00,
 -5.06914062e+05,  0.00000000e+00,  0.00000000e+00,
 -5.06914062e+05,  0.00000000e+00,  0.00000000e+00,
  5.45154948e+05,  0.00000000e+00,  0.00000000e+00,
  5.19364583e+05,  0.00000000e+00,  0.00000000e+00,
  1.50829167e+06,  0.00000000e+00,  0.00000000e+00,
  1.50829167e+06,  0.00000000e+00,  0.00000000e+00,
  1.80888281e+06,  0.00000000e+00,  0.00000000e+00,
  1.12054688e+06,  0.00000000e+00,  0.00000000e+00,
  1.12054688e+06,  0.00000000e+00,  0.00000000e+00,
  0.00000000e+00,  1.24238411e+06,  0.00000000e+00,
  0.00000000e+00, -5.06914062e+05,  0.00000000e+00,
  0.00000000e+00, -5.06914062e+05,  0.00000000e+00,
  0.00000000e+00,  5.45154948e+05,  0.00000000e+00,
  0.00000000e+00,  5.19364583e+05,  0.00000000e+00,
  0.00000000e+00,  1.50829167e+06,  0.00000000e+00,
  0.00000000e+00,  1.50829167e+06,  0.00000000e+00,
  0.00000000e+00,  1.80888281e+06,  0.00000000e+00,
  0.00000000e+00,  1.12054688e+06,  0.00000000e+00,
  0.00000000e+00,  1.12054688e+06,  0.00000000e+00,
  0.00000000e+00,  0.00000000e+00,  1.24238411e+06,
  0.00000000e+00,  0.00000000e+00, -5.06914062e+05,
  0.00000000e+00,  0.00000000e+00, -5.06914062e+05,
  0.00000000e+00,  0.00000000e+00,  5.45154948e+05,
  0.00000000e+00,  0.00000000e+00,  5.19364583e+05,
  0.00000000e+00,  0.00000000e+00,  1.50829167e+06,
  0.00000000e+00,  0.00000000e+00,  1.50829167e+06,
  0.00000000e+00,  0.00000000e+00,  1.80888281e+06,
  0.00000000e+00,  0.00000000e+00,  1.12054688e+06,
  0.00000000e+00,  0.00000000e+00,  1.12054688e+06,
 -3.13041667e+05,  0.00000000e+00,  0.00000000e+00,
  1.46026823e+06,  0.00000000e+00,  0.00000000e+00,
 -6.43869792e+05,  0.00000000e+00,  0.00000000e+00,
  1.46026823e+06,  0.00000000e+00,  0.00000000e+00,
  2.05611458e+06,  0.00000000e+00,  0.00000000e+00,
  5.19364583e+05,  0.00000000e+00,  0.00000000e+00,
  1.12054688e+06,  0.00000000e+00,  0.00000000e+00,
  1.12054688e+06,  0.00000000e+00,  0.00000000e+00,
  3.25847917e+06,  0.00000000e+00,  0.00000000e+00,
  2.05611458e+06,  0.00000000e+00,  0.00000000e+00,
  0.00000000e+00, -3.13041667e+05,  0.00000000e+00,
  0.00000000e+00,  1.46026823e+06,  0.00000000e+00,
  0.00000000e+00, -6.43869792e+05,  0.00000000e+00,
  0.00000000e+00,  1.46026823e+06,  0.00000000e+00,
  0.00000000e+00,  2.05611458e+06,  0.00000000e+00,
  0.00000000e+00,  5.19364583e+05,  0.00000000e+00,
  0.00000000e+00,  1.12054688e+06,  0.00000000e+00,
  0.00000000e+00,  1.12054688e+06,  0.00000000e+00,
  0.00000000e+00,  3.25847917e+06,  0.00000000e+00,
  0.00000000e+00,  2.05611458e+06,  0.00000000e+00,
  0.00000000e+00,  0.00000000e+00, -3.13041667e+05,
  0.00000000e+00,  0.00000000e+00,  1.46026823e+06,
  0.00000000e+00,  0.00000000e+00, -6.43869792e+05,
  0.00000000e+00,  0.00000000e+00,  1.46026823e+06,
  0.00000000e+00,  0.00000000e+00,  2.05611458e+06,
  0.00000000e+00,  0.00000000e+00,  5.19364583e+05,
  0.00000000e+00,  0.00000000e+00,  1.12054688e+06,
  0.00000000e+00,  0.00000000e+00,  1.12054688e+06,
  0.00000000e+00,  0.00000000e+00,  3.25847917e+06,
  0.00000000e+00,  0.00000000e+00,  2.05611458e+06,
 -3.13041667e+05,  0.00000000e+00,  0.00000000e+00,
 -6.43869792e+05,  0.00000000e+00,  0.00000000e+00,
  1.46026823e+06,  0.00000000e+00,  0.00000000e+00,
  1.46026823e+06,  0.00000000e+00,  0.00000000e+00,
  2.05611458e+06,  0.00000000e+00,  0.00000000e+00,
  1.12054688e+06,  0.00000000e+00,  0.00000000e+00,
  5.19364583e+05,  0.00000000e+00,  0.00000000e+00,
  1.12054688e+06,  0.00000000e+00,  0.00000000e+00,
  2.05611458e+06,  0.00000000e+00,  0.00000000e+00,
  3.25847917e+06,  0.00000000e+00,  0.00000000e+00,
  0.00000000e+00, -3.13041667e+05,  0.00000000e+00,
  0.00000000e+00, -6.43869792e+05,  0.00000000e+00,
  0.00000000e+00,  1.46026823e+06,  0.00000000e+00,
  0.00000000e+00,  1.46026823e+06,  0.00000000e+00,
  0.00000000e+00,  2.05611458e+06,  0.00000000e+00,
  0.00000000e+00,  1.12054688e+06,  0.00000000e+00,
  0.00000000e+00,  5.19364583e+05,  0.00000000e+00,
  0.00000000e+00,  1.12054688e+06,  0.00000000e+00,
  0.00000000e+00,  2.05611458e+06,  0.00000000e+00,
  0.00000000e+00,  3.25847917e+06,  0.00000000e+00,
  0.00000000e+00,  0.00000000e+00, -3.13041667e+05,
  0.00000000e+00,  0.00000000e+00, -6.43869792e+05,
  0.00000000e+00,  0.00000000e+00,  1.46026823e+06,
  0.00000000e+00,  0.00000000e+00,  1.46026823e+06,
  0.00000000e+00,  0.00000000e+00,  2.05611458e+06,
  0.00000000e+00,  0.00000000e+00,  1.12054688e+06,
  0.00000000e+00,  0.00000000e+00,  5.19364583e+05,
  0.00000000e+00,  0.00000000e+00,  1.12054688e+06,
  0.00000000e+00,  0.00000000e+00,  2.05611458e+06,
  0.00000000e+00,  0.00000000e+00,  3.25847917e+06,
};

const double pylith::feassemble::ElasticityExplicitGravData3DQuadratic::_valsResidualLumped[] = {
  2.17441536e+10, -9.08020975e+09, -5.09679086e+09,
  6.22351054e+10, -2.47908763e+10, -1.65890189e+10,
 -4.60097182e+10, -5.34117225e+10, -5.55302502e+10,
 -1.00516830e+10,  4.19873726e+10, -6.27434896e+10,
 -8.42078211e+09,  6.33182361e+10, -1.24376360e+11,
  6.62278190e+10,  1.09244311e+11, -1.20869461e+11,
 -9.59516349e+10, -5.29561314e+10, -1.08713086e+11,
 -1.25237479e+10, -6.06732213e+10, -7.32158278e+10,
 -3.24416221e+10,  5.14211527e+10, -1.13633369e+11,
  5.52060281e+10, -6.50368653e+10, -3.06762643e+10,
};

const double pylith::feassemble::ElasticityExplicitGravData3DQuadratic::_valsJacobianLumped[] = {
  2.57903646e+06,  2.57903646e+06,  2.57903646e+06,
  2.40117188e+06,  2.40117188e+06,  2.40117188e+06,
  2.40117188e+06,  2.40117188e+06,  2.40117188e+06,
  2.40117188e+06,  2.40117188e+06,  2.40117188e+06,
  1.20947917e+07,  1.20947917e+07,  1.20947917e+07,
  8.35963542e+06,  8.35963542e+06,  8.35963542e+06,
  8.35963542e+06,  8.35963542e+06,  8.35963542e+06,
  8.35963542e+06,  8.35963542e+06,  8.35963542e+06,
  1.20947917e+07,  1.20947917e+07,  1.20947917e+07,
  1.20947917e+07,  1.20947917e+07,  1.20947917e+07,
};

pylith::feassemble::ElasticityExplicitGravData3DQuadratic::ElasticityExplicitGravData3DQuadratic(void)
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
  fieldTIncr = const_cast<double*>(_fieldTIncr);
  fieldT = const_cast<double*>(_fieldT);
  fieldTmdt = const_cast<double*>(_fieldTmdt);
  valsResidual = const_cast<double*>(_valsResidual);
  valsJacobian = const_cast<double*>(_valsJacobian);
  valsResidualLumped = const_cast<double*>(_valsResidualLumped);
  valsJacobianLumped = const_cast<double*>(_valsJacobianLumped);
} // constructor

pylith::feassemble::ElasticityExplicitGravData3DQuadratic::~ElasticityExplicitGravData3DQuadratic(void)
{}


// End of file
