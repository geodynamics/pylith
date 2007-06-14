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

#include "ElasticityExplicitData3DQuadratic.hh"

const int pylith::feassemble::ElasticityExplicitData3DQuadratic::_spaceDim = 3;

const int pylith::feassemble::ElasticityExplicitData3DQuadratic::_cellDim = 3;

const int pylith::feassemble::ElasticityExplicitData3DQuadratic::_numVertices = 10;

const int pylith::feassemble::ElasticityExplicitData3DQuadratic::_numCells = 1;

const int pylith::feassemble::ElasticityExplicitData3DQuadratic::_numBasis = 10;

const int pylith::feassemble::ElasticityExplicitData3DQuadratic::_numQuadPts = 4;

const char* pylith::feassemble::ElasticityExplicitData3DQuadratic::_matType = "ElasticIsotropic3D";

const char* pylith::feassemble::ElasticityExplicitData3DQuadratic::_matDBFilename = "data/elasticisotropic3d.spatialdb";

const int pylith::feassemble::ElasticityExplicitData3DQuadratic::_matId = 0;

const char* pylith::feassemble::ElasticityExplicitData3DQuadratic::_matLabel = "elastic isotropic 3-D";

const double pylith::feassemble::ElasticityExplicitData3DQuadratic::_dt =   1.00000000e-02;

const double pylith::feassemble::ElasticityExplicitData3DQuadratic::_vertices[] = {
 -5.00000000e-01, -2.00000000e+00, -1.00000000e+00,
  2.00000000e+00, -2.00000000e+00, -5.00000000e-01,
  1.00000000e+00,  1.00000000e+00,  0.00000000e+00,
  2.00000000e-01,  5.00000000e-01,  2.00000000e+00,
  7.00000000e-01, -2.10000000e+00, -8.00000000e-01,
  3.00000000e-01, -5.00000000e-01, -5.00000000e-01,
 -2.00000000e-01, -8.00000000e-01,  5.00000000e-01,
  1.50000000e+00, -6.00000000e-01, -2.00000000e-01,
  6.00000000e-01,  8.00000000e-01,  9.00000000e-01,
  1.10000000e+00, -8.00000000e-01,  7.00000000e-01,
};

const int pylith::feassemble::ElasticityExplicitData3DQuadratic::_cells[] = {
0,1,2,3,4,5,6,7,8,9,
};

const double pylith::feassemble::ElasticityExplicitData3DQuadratic::_verticesRef[] = {
 -1.00000000e+00, -1.00000000e+00, -1.00000000e+00,
  1.00000000e+00, -1.00000000e+00, -1.00000000e+00,
 -1.00000000e+00,  1.00000000e+00, -1.00000000e+00,
 -1.00000000e+00, -1.00000000e+00,  1.00000000e+00,
 -1.00000000e+00,  0.00000000e+00, -1.00000000e+00,
  0.00000000e+00,  0.00000000e+00, -1.00000000e+00,
  0.00000000e+00, -1.00000000e+00, -1.00000000e+00,
 -1.00000000e+00, -1.00000000e+00,  0.00000000e+00,
  0.00000000e+00, -1.00000000e+00,  0.00000000e+00,
 -1.00000000e+00,  0.00000000e+00,  0.00000000e+00,
};

const double pylith::feassemble::ElasticityExplicitData3DQuadratic::_quadPts[] = {
  8.33333333e-02,  8.33333333e-02,  8.33333333e-02,
  7.50000000e-01,  8.33333333e-02,  8.33333333e-02,
  8.33333333e-02,  7.50000000e-01,  8.33333333e-02,
  8.33333333e-02,  8.33333333e-02,  7.50000000e-01,
};

const double pylith::feassemble::ElasticityExplicitData3DQuadratic::_quadWts[] = {
  1.25000000e-01,  1.25000000e-01,  1.25000000e-01,  1.25000000e-01,
};

const double pylith::feassemble::ElasticityExplicitData3DQuadratic::_basis[] = {
  3.75000000e-01, -6.94444444e-02, -6.94444444e-02,
 -6.94444444e-02,  2.50000000e-01,  2.50000000e-01,
  2.50000000e-01,  2.77777778e-02,  2.77777778e-02,
  2.77777778e-02, -6.94444444e-02,  3.75000000e-01,
 -6.94444444e-02, -6.94444444e-02,  2.50000000e-01,
  2.77777778e-02,  2.77777778e-02,  2.50000000e-01,
  2.77777778e-02,  2.50000000e-01, -6.94444444e-02,
 -6.94444444e-02,  3.75000000e-01, -6.94444444e-02,
  2.77777778e-02,  2.50000000e-01,  2.77777778e-02,
  2.50000000e-01,  2.50000000e-01,  2.77777778e-02,
 -6.94444444e-02, -6.94444444e-02, -6.94444444e-02,
  3.75000000e-01,  2.77777778e-02,  2.77777778e-02,
  2.50000000e-01,  2.77777778e-02,  2.50000000e-01,
  2.50000000e-01,};

const double pylith::feassemble::ElasticityExplicitData3DQuadratic::_basisDerivRef[] = {
 -2.00000000e+00, -2.00000000e+00, -2.00000000e+00,
 -6.66666667e-01,  0.00000000e+00,  0.00000000e+00,
  0.00000000e+00, -6.66666667e-01,  0.00000000e+00,
  0.00000000e+00,  0.00000000e+00, -6.66666667e-01,
  2.66666667e+00, -3.33333333e-01, -3.33333333e-01,
 -3.33333333e-01,  2.66666667e+00, -3.33333333e-01,
 -3.33333333e-01, -3.33333333e-01,  2.66666667e+00,
  3.33333333e-01,  3.33333333e-01,  0.00000000e+00,
  0.00000000e+00,  3.33333333e-01,  3.33333333e-01,
  3.33333333e-01,  0.00000000e+00,  3.33333333e-01,
  6.66666667e-01,  6.66666667e-01,  6.66666667e-01,
  2.00000000e+00,  0.00000000e+00,  0.00000000e+00,
  0.00000000e+00, -6.66666667e-01,  0.00000000e+00,
  0.00000000e+00,  0.00000000e+00, -6.66666667e-01,
 -2.66666667e+00, -3.00000000e+00, -3.00000000e+00,
 -3.33333333e-01,  1.11022302e-16, -3.33333333e-01,
 -3.33333333e-01, -3.33333333e-01,  1.11022302e-16,
  3.33333333e-01,  3.00000000e+00,  0.00000000e+00,
  0.00000000e+00,  3.33333333e-01,  3.33333333e-01,
  3.33333333e-01,  0.00000000e+00,  3.00000000e+00,
  6.66666667e-01,  6.66666667e-01,  6.66666667e-01,
 -6.66666667e-01,  0.00000000e+00,  0.00000000e+00,
  0.00000000e+00,  2.00000000e+00,  0.00000000e+00,
  0.00000000e+00,  0.00000000e+00, -6.66666667e-01,
  1.11022302e-16, -3.33333333e-01, -3.33333333e-01,
 -3.00000000e+00, -2.66666667e+00, -3.00000000e+00,
 -3.33333333e-01, -3.33333333e-01, -1.11022302e-16,
  3.00000000e+00,  3.33333333e-01,  0.00000000e+00,
  0.00000000e+00,  3.33333333e-01,  3.00000000e+00,
  3.33333333e-01,  0.00000000e+00,  3.33333333e-01,
  6.66666667e-01,  6.66666667e-01,  6.66666667e-01,
 -6.66666667e-01,  0.00000000e+00,  0.00000000e+00,
  0.00000000e+00, -6.66666667e-01,  0.00000000e+00,
  0.00000000e+00,  0.00000000e+00,  2.00000000e+00,
 -1.11022302e-16, -3.33333333e-01, -3.33333333e-01,
 -3.33333333e-01, -1.11022302e-16, -3.33333333e-01,
 -3.00000000e+00, -3.00000000e+00, -2.66666667e+00,
  3.33333333e-01,  3.33333333e-01,  0.00000000e+00,
  0.00000000e+00,  3.00000000e+00,  3.33333333e-01,
  3.00000000e+00,  0.00000000e+00,  3.33333333e-01,
};

const double pylith::feassemble::ElasticityExplicitData3DQuadratic::_fieldTpdt[] = {
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

const double pylith::feassemble::ElasticityExplicitData3DQuadratic::_fieldT[] = {
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

const double pylith::feassemble::ElasticityExplicitData3DQuadratic::_fieldTmdt[] = {
  2.00000000e-01, -3.00000000e-01, -4.00000000e-01,
 -6.00000000e-01,  2.00000000e-01,  3.00000000e-01,
  5.00000000e-01,  2.00000000e-01,  5.00000000e-01,
 -3.00000000e-01, -6.00000000e-01, -3.00000000e-01,
 -5.00000000e-01, -9.00000000e-01,  4.00000000e-01,
 -3.00000000e-01, -6.00000000e-01, -8.00000000e-01,
  9.00000000e-01,  5.00000000e-01, -2.00000000e-01,
 -7.00000000e-01, -3.00000000e-01, -3.00000000e-01,
 -5.00000000e-01, -8.00000000e-01,  4.00000000e-01,
 -4.00000000e-01,  5.00000000e-01, -8.00000000e-01,
};

const double pylith::feassemble::ElasticityExplicitData3DQuadratic::_valsResidual[] = {
  8.00450011e+10, -9.00472985e+09,  1.61541766e+11,
 -3.54088543e+10, -6.47782824e+10,  2.41769691e+10,
 -1.55663344e+11, -4.18485549e+11,  2.08248421e+11,
 -9.98985902e+09, -4.26541080e+10,  3.19516877e+11,
 -5.76612056e+10, -2.26094787e+11, -1.40838144e+11,
  1.37516193e+11,  4.40486089e+11, -1.93421106e+11,
 -8.06475029e+11, -3.08666160e+11, -1.91097837e+11,
  8.66282546e+10,  4.68394314e+11, -2.59144518e+11,
  1.66627489e+11,  6.19155802e+11, -4.59691032e+11,
  5.94347823e+11, -4.58318391e+11,  5.30727623e+11,
};

const double pylith::feassemble::ElasticityExplicitData3DQuadratic::_valsJacobian[] = {
  8.58871421e+06,  0.00000000e+00,  0.00000000e+00,
 -2.12223687e+06,  0.00000000e+00,  0.00000000e+00,
 -2.14181277e+06,  0.00000000e+00,  0.00000000e+00,
 -2.72868977e+06,  0.00000000e+00,  0.00000000e+00,
  4.16422397e+06,  0.00000000e+00,  0.00000000e+00,
  4.15443601e+06,  0.00000000e+00,  0.00000000e+00,
  3.86099751e+06,  0.00000000e+00,  0.00000000e+00,
 -1.20103952e+06,  0.00000000e+00,  0.00000000e+00,
 -1.50426598e+06,  0.00000000e+00,  0.00000000e+00,
 -1.49447802e+06,  0.00000000e+00,  0.00000000e+00,
  0.00000000e+00,  8.58871421e+06,  0.00000000e+00,
  0.00000000e+00, -2.12223687e+06,  0.00000000e+00,
  0.00000000e+00, -2.14181277e+06,  0.00000000e+00,
  0.00000000e+00, -2.72868977e+06,  0.00000000e+00,
  0.00000000e+00,  4.16422397e+06,  0.00000000e+00,
  0.00000000e+00,  4.15443601e+06,  0.00000000e+00,
  0.00000000e+00,  3.86099751e+06,  0.00000000e+00,
  0.00000000e+00, -1.20103952e+06,  0.00000000e+00,
  0.00000000e+00, -1.50426598e+06,  0.00000000e+00,
  0.00000000e+00, -1.49447802e+06,  0.00000000e+00,
  0.00000000e+00,  0.00000000e+00,  8.58871421e+06,
  0.00000000e+00,  0.00000000e+00, -2.12223687e+06,
  0.00000000e+00,  0.00000000e+00, -2.14181277e+06,
  0.00000000e+00,  0.00000000e+00, -2.72868977e+06,
  0.00000000e+00,  0.00000000e+00,  4.16422397e+06,
  0.00000000e+00,  0.00000000e+00,  4.15443601e+06,
  0.00000000e+00,  0.00000000e+00,  3.86099751e+06,
  0.00000000e+00,  0.00000000e+00, -1.20103952e+06,
  0.00000000e+00,  0.00000000e+00, -1.50426598e+06,
  0.00000000e+00,  0.00000000e+00, -1.49447802e+06,
 -2.12223687e+06,  0.00000000e+00,  0.00000000e+00,
  7.46054777e+06,  0.00000000e+00,  0.00000000e+00,
 -1.88541131e+06,  0.00000000e+00,  0.00000000e+00,
 -2.47228831e+06,  0.00000000e+00,  0.00000000e+00,
  3.24117870e+06,  0.00000000e+00,  0.00000000e+00,
 -1.43180084e+06,  0.00000000e+00,  0.00000000e+00,
 -1.72523934e+06,  0.00000000e+00,  0.00000000e+00,
  3.35959148e+06,  0.00000000e+00,  0.00000000e+00,
 -1.60682656e+06,  0.00000000e+00,  0.00000000e+00,
  3.06615298e+06,  0.00000000e+00,  0.00000000e+00,
  0.00000000e+00, -2.12223687e+06,  0.00000000e+00,
  0.00000000e+00,  7.46054777e+06,  0.00000000e+00,
  0.00000000e+00, -1.88541131e+06,  0.00000000e+00,
  0.00000000e+00, -2.47228831e+06,  0.00000000e+00,
  0.00000000e+00,  3.24117870e+06,  0.00000000e+00,
  0.00000000e+00, -1.43180084e+06,  0.00000000e+00,
  0.00000000e+00, -1.72523934e+06,  0.00000000e+00,
  0.00000000e+00,  3.35959148e+06,  0.00000000e+00,
  0.00000000e+00, -1.60682656e+06,  0.00000000e+00,
  0.00000000e+00,  3.06615298e+06,  0.00000000e+00,
  0.00000000e+00,  0.00000000e+00, -2.12223687e+06,
  0.00000000e+00,  0.00000000e+00,  7.46054777e+06,
  0.00000000e+00,  0.00000000e+00, -1.88541131e+06,
  0.00000000e+00,  0.00000000e+00, -2.47228831e+06,
  0.00000000e+00,  0.00000000e+00,  3.24117870e+06,
  0.00000000e+00,  0.00000000e+00, -1.43180084e+06,
  0.00000000e+00,  0.00000000e+00, -1.72523934e+06,
  0.00000000e+00,  0.00000000e+00,  3.35959148e+06,
  0.00000000e+00,  0.00000000e+00, -1.60682656e+06,
  0.00000000e+00,  0.00000000e+00,  3.06615298e+06,
 -2.14181277e+06,  0.00000000e+00,  0.00000000e+00,
 -1.88541131e+06,  0.00000000e+00,  0.00000000e+00,
  7.54668174e+06,  0.00000000e+00,  0.00000000e+00,
 -2.49186421e+06,  0.00000000e+00,  0.00000000e+00,
 -1.41418253e+06,  0.00000000e+00,  0.00000000e+00,
  3.30186400e+06,  0.00000000e+00,  0.00000000e+00,
 -1.71740898e+06,  0.00000000e+00,  0.00000000e+00,
  3.43006473e+06,  0.00000000e+00,  0.00000000e+00,
  3.12683828e+06,  0.00000000e+00,  0.00000000e+00,
 -1.58920825e+06,  0.00000000e+00,  0.00000000e+00,
  0.00000000e+00, -2.14181277e+06,  0.00000000e+00,
  0.00000000e+00, -1.88541131e+06,  0.00000000e+00,
  0.00000000e+00,  7.54668174e+06,  0.00000000e+00,
  0.00000000e+00, -2.49186421e+06,  0.00000000e+00,
  0.00000000e+00, -1.41418253e+06,  0.00000000e+00,
  0.00000000e+00,  3.30186400e+06,  0.00000000e+00,
  0.00000000e+00, -1.71740898e+06,  0.00000000e+00,
  0.00000000e+00,  3.43006473e+06,  0.00000000e+00,
  0.00000000e+00,  3.12683828e+06,  0.00000000e+00,
  0.00000000e+00, -1.58920825e+06,  0.00000000e+00,
  0.00000000e+00,  0.00000000e+00, -2.14181277e+06,
  0.00000000e+00,  0.00000000e+00, -1.88541131e+06,
  0.00000000e+00,  0.00000000e+00,  7.54668174e+06,
  0.00000000e+00,  0.00000000e+00, -2.49186421e+06,
  0.00000000e+00,  0.00000000e+00, -1.41418253e+06,
  0.00000000e+00,  0.00000000e+00,  3.30186400e+06,
  0.00000000e+00,  0.00000000e+00, -1.71740898e+06,
  0.00000000e+00,  0.00000000e+00,  3.43006473e+06,
  0.00000000e+00,  0.00000000e+00,  3.12683828e+06,
  0.00000000e+00,  0.00000000e+00, -1.58920825e+06,
 -2.72868977e+06,  0.00000000e+00,  0.00000000e+00,
 -2.47228831e+06,  0.00000000e+00,  0.00000000e+00,
 -2.49186421e+06,  0.00000000e+00,  0.00000000e+00,
  1.01289405e+07,  0.00000000e+00,  0.00000000e+00,
 -1.17943173e+06,  0.00000000e+00,  0.00000000e+00,
 -1.18921968e+06,  0.00000000e+00,  0.00000000e+00,
  5.12118270e+06,  0.00000000e+00,  0.00000000e+00,
 -1.06101895e+06,  0.00000000e+00,  0.00000000e+00,
  5.23959548e+06,  0.00000000e+00,  0.00000000e+00,
  5.24938343e+06,  0.00000000e+00,  0.00000000e+00,
  0.00000000e+00, -2.72868977e+06,  0.00000000e+00,
  0.00000000e+00, -2.47228831e+06,  0.00000000e+00,
  0.00000000e+00, -2.49186421e+06,  0.00000000e+00,
  0.00000000e+00,  1.01289405e+07,  0.00000000e+00,
  0.00000000e+00, -1.17943173e+06,  0.00000000e+00,
  0.00000000e+00, -1.18921968e+06,  0.00000000e+00,
  0.00000000e+00,  5.12118270e+06,  0.00000000e+00,
  0.00000000e+00, -1.06101895e+06,  0.00000000e+00,
  0.00000000e+00,  5.23959548e+06,  0.00000000e+00,
  0.00000000e+00,  5.24938343e+06,  0.00000000e+00,
  0.00000000e+00,  0.00000000e+00, -2.72868977e+06,
  0.00000000e+00,  0.00000000e+00, -2.47228831e+06,
  0.00000000e+00,  0.00000000e+00, -2.49186421e+06,
  0.00000000e+00,  0.00000000e+00,  1.01289405e+07,
  0.00000000e+00,  0.00000000e+00, -1.17943173e+06,
  0.00000000e+00,  0.00000000e+00, -1.18921968e+06,
  0.00000000e+00,  0.00000000e+00,  5.12118270e+06,
  0.00000000e+00,  0.00000000e+00, -1.06101895e+06,
  0.00000000e+00,  0.00000000e+00,  5.23959548e+06,
  0.00000000e+00,  0.00000000e+00,  5.24938343e+06,
  4.16422397e+06,  0.00000000e+00,  0.00000000e+00,
  3.24117870e+06,  0.00000000e+00,  0.00000000e+00,
 -1.41418253e+06,  0.00000000e+00,  0.00000000e+00,
 -1.17943173e+06,  0.00000000e+00,  0.00000000e+00,
  6.50957790e+06,  0.00000000e+00,  0.00000000e+00,
  4.18189729e+06,  0.00000000e+00,  0.00000000e+00,
  4.29927269e+06,  0.00000000e+00,  0.00000000e+00,
  3.72037466e+06,  0.00000000e+00,  0.00000000e+00,
  1.51006944e+06,  0.00000000e+00,  0.00000000e+00,
  3.83775006e+06,  0.00000000e+00,  0.00000000e+00,
  0.00000000e+00,  4.16422397e+06,  0.00000000e+00,
  0.00000000e+00,  3.24117870e+06,  0.00000000e+00,
  0.00000000e+00, -1.41418253e+06,  0.00000000e+00,
  0.00000000e+00, -1.17943173e+06,  0.00000000e+00,
  0.00000000e+00,  6.50957790e+06,  0.00000000e+00,
  0.00000000e+00,  4.18189729e+06,  0.00000000e+00,
  0.00000000e+00,  4.29927269e+06,  0.00000000e+00,
  0.00000000e+00,  3.72037466e+06,  0.00000000e+00,
  0.00000000e+00,  1.51006944e+06,  0.00000000e+00,
  0.00000000e+00,  3.83775006e+06,  0.00000000e+00,
  0.00000000e+00,  0.00000000e+00,  4.16422397e+06,
  0.00000000e+00,  0.00000000e+00,  3.24117870e+06,
  0.00000000e+00,  0.00000000e+00, -1.41418253e+06,
  0.00000000e+00,  0.00000000e+00, -1.17943173e+06,
  0.00000000e+00,  0.00000000e+00,  6.50957790e+06,
  0.00000000e+00,  0.00000000e+00,  4.18189729e+06,
  0.00000000e+00,  0.00000000e+00,  4.29927269e+06,
  0.00000000e+00,  0.00000000e+00,  3.72037466e+06,
  0.00000000e+00,  0.00000000e+00,  1.51006944e+06,
  0.00000000e+00,  0.00000000e+00,  3.83775006e+06,
  4.15443601e+06,  0.00000000e+00,  0.00000000e+00,
 -1.43180084e+06,  0.00000000e+00,  0.00000000e+00,
  3.30186400e+06,  0.00000000e+00,  0.00000000e+00,
 -1.18921968e+06,  0.00000000e+00,  0.00000000e+00,
  4.18189729e+06,  0.00000000e+00,  0.00000000e+00,
  6.54872971e+06,  0.00000000e+00,  0.00000000e+00,
  4.30318787e+06,  0.00000000e+00,  0.00000000e+00,
  3.75561128e+06,  0.00000000e+00,  0.00000000e+00,
  3.87690186e+06,  0.00000000e+00,  0.00000000e+00,
  1.51006944e+06,  0.00000000e+00,  0.00000000e+00,
  0.00000000e+00,  4.15443601e+06,  0.00000000e+00,
  0.00000000e+00, -1.43180084e+06,  0.00000000e+00,
  0.00000000e+00,  3.30186400e+06,  0.00000000e+00,
  0.00000000e+00, -1.18921968e+06,  0.00000000e+00,
  0.00000000e+00,  4.18189729e+06,  0.00000000e+00,
  0.00000000e+00,  6.54872971e+06,  0.00000000e+00,
  0.00000000e+00,  4.30318787e+06,  0.00000000e+00,
  0.00000000e+00,  3.75561128e+06,  0.00000000e+00,
  0.00000000e+00,  3.87690186e+06,  0.00000000e+00,
  0.00000000e+00,  1.51006944e+06,  0.00000000e+00,
  0.00000000e+00,  0.00000000e+00,  4.15443601e+06,
  0.00000000e+00,  0.00000000e+00, -1.43180084e+06,
  0.00000000e+00,  0.00000000e+00,  3.30186400e+06,
  0.00000000e+00,  0.00000000e+00, -1.18921968e+06,
  0.00000000e+00,  0.00000000e+00,  4.18189729e+06,
  0.00000000e+00,  0.00000000e+00,  6.54872971e+06,
  0.00000000e+00,  0.00000000e+00,  4.30318787e+06,
  0.00000000e+00,  0.00000000e+00,  3.75561128e+06,
  0.00000000e+00,  0.00000000e+00,  3.87690186e+06,
  0.00000000e+00,  0.00000000e+00,  1.51006944e+06,
  3.86099751e+06,  0.00000000e+00,  0.00000000e+00,
 -1.72523934e+06,  0.00000000e+00,  0.00000000e+00,
 -1.71740898e+06,  0.00000000e+00,  0.00000000e+00,
  5.12118270e+06,  0.00000000e+00,  0.00000000e+00,
  4.29927269e+06,  0.00000000e+00,  0.00000000e+00,
  4.30318787e+06,  0.00000000e+00,  0.00000000e+00,
  7.72248371e+06,  0.00000000e+00,  0.00000000e+00,
  1.51006944e+06,  0.00000000e+00,  0.00000000e+00,
  4.93328046e+06,  0.00000000e+00,  0.00000000e+00,
  4.92936528e+06,  0.00000000e+00,  0.00000000e+00,
  0.00000000e+00,  3.86099751e+06,  0.00000000e+00,
  0.00000000e+00, -1.72523934e+06,  0.00000000e+00,
  0.00000000e+00, -1.71740898e+06,  0.00000000e+00,
  0.00000000e+00,  5.12118270e+06,  0.00000000e+00,
  0.00000000e+00,  4.29927269e+06,  0.00000000e+00,
  0.00000000e+00,  4.30318787e+06,  0.00000000e+00,
  0.00000000e+00,  7.72248371e+06,  0.00000000e+00,
  0.00000000e+00,  1.51006944e+06,  0.00000000e+00,
  0.00000000e+00,  4.93328046e+06,  0.00000000e+00,
  0.00000000e+00,  4.92936528e+06,  0.00000000e+00,
  0.00000000e+00,  0.00000000e+00,  3.86099751e+06,
  0.00000000e+00,  0.00000000e+00, -1.72523934e+06,
  0.00000000e+00,  0.00000000e+00, -1.71740898e+06,
  0.00000000e+00,  0.00000000e+00,  5.12118270e+06,
  0.00000000e+00,  0.00000000e+00,  4.29927269e+06,
  0.00000000e+00,  0.00000000e+00,  4.30318787e+06,
  0.00000000e+00,  0.00000000e+00,  7.72248371e+06,
  0.00000000e+00,  0.00000000e+00,  1.51006944e+06,
  0.00000000e+00,  0.00000000e+00,  4.93328046e+06,
  0.00000000e+00,  0.00000000e+00,  4.92936528e+06,
 -1.20103952e+06,  0.00000000e+00,  0.00000000e+00,
  3.35959148e+06,  0.00000000e+00,  0.00000000e+00,
  3.43006473e+06,  0.00000000e+00,  0.00000000e+00,
 -1.06101895e+06,  0.00000000e+00,  0.00000000e+00,
  3.72037466e+06,  0.00000000e+00,  0.00000000e+00,
  3.75561128e+06,  0.00000000e+00,  0.00000000e+00,
  1.51006944e+06,  0.00000000e+00,  0.00000000e+00,
  6.03592678e+06,  0.00000000e+00,  0.00000000e+00,
  3.82562157e+06,  0.00000000e+00,  0.00000000e+00,
  3.79038495e+06,  0.00000000e+00,  0.00000000e+00,
  0.00000000e+00, -1.20103952e+06,  0.00000000e+00,
  0.00000000e+00,  3.35959148e+06,  0.00000000e+00,
  0.00000000e+00,  3.43006473e+06,  0.00000000e+00,
  0.00000000e+00, -1.06101895e+06,  0.00000000e+00,
  0.00000000e+00,  3.72037466e+06,  0.00000000e+00,
  0.00000000e+00,  3.75561128e+06,  0.00000000e+00,
  0.00000000e+00,  1.51006944e+06,  0.00000000e+00,
  0.00000000e+00,  6.03592678e+06,  0.00000000e+00,
  0.00000000e+00,  3.82562157e+06,  0.00000000e+00,
  0.00000000e+00,  3.79038495e+06,  0.00000000e+00,
  0.00000000e+00,  0.00000000e+00, -1.20103952e+06,
  0.00000000e+00,  0.00000000e+00,  3.35959148e+06,
  0.00000000e+00,  0.00000000e+00,  3.43006473e+06,
  0.00000000e+00,  0.00000000e+00, -1.06101895e+06,
  0.00000000e+00,  0.00000000e+00,  3.72037466e+06,
  0.00000000e+00,  0.00000000e+00,  3.75561128e+06,
  0.00000000e+00,  0.00000000e+00,  1.51006944e+06,
  0.00000000e+00,  0.00000000e+00,  6.03592678e+06,
  0.00000000e+00,  0.00000000e+00,  3.82562157e+06,
  0.00000000e+00,  0.00000000e+00,  3.79038495e+06,
 -1.50426598e+06,  0.00000000e+00,  0.00000000e+00,
 -1.60682656e+06,  0.00000000e+00,  0.00000000e+00,
  3.12683828e+06,  0.00000000e+00,  0.00000000e+00,
  5.23959548e+06,  0.00000000e+00,  0.00000000e+00,
  1.51006944e+06,  0.00000000e+00,  0.00000000e+00,
  3.87690186e+06,  0.00000000e+00,  0.00000000e+00,
  4.93328046e+06,  0.00000000e+00,  0.00000000e+00,
  3.82562157e+06,  0.00000000e+00,  0.00000000e+00,
  7.24883259e+06,  0.00000000e+00,  0.00000000e+00,
  4.88200017e+06,  0.00000000e+00,  0.00000000e+00,
  0.00000000e+00, -1.50426598e+06,  0.00000000e+00,
  0.00000000e+00, -1.60682656e+06,  0.00000000e+00,
  0.00000000e+00,  3.12683828e+06,  0.00000000e+00,
  0.00000000e+00,  5.23959548e+06,  0.00000000e+00,
  0.00000000e+00,  1.51006944e+06,  0.00000000e+00,
  0.00000000e+00,  3.87690186e+06,  0.00000000e+00,
  0.00000000e+00,  4.93328046e+06,  0.00000000e+00,
  0.00000000e+00,  3.82562157e+06,  0.00000000e+00,
  0.00000000e+00,  7.24883259e+06,  0.00000000e+00,
  0.00000000e+00,  4.88200017e+06,  0.00000000e+00,
  0.00000000e+00,  0.00000000e+00, -1.50426598e+06,
  0.00000000e+00,  0.00000000e+00, -1.60682656e+06,
  0.00000000e+00,  0.00000000e+00,  3.12683828e+06,
  0.00000000e+00,  0.00000000e+00,  5.23959548e+06,
  0.00000000e+00,  0.00000000e+00,  1.51006944e+06,
  0.00000000e+00,  0.00000000e+00,  3.87690186e+06,
  0.00000000e+00,  0.00000000e+00,  4.93328046e+06,
  0.00000000e+00,  0.00000000e+00,  3.82562157e+06,
  0.00000000e+00,  0.00000000e+00,  7.24883259e+06,
  0.00000000e+00,  0.00000000e+00,  4.88200017e+06,
 -1.49447802e+06,  0.00000000e+00,  0.00000000e+00,
  3.06615298e+06,  0.00000000e+00,  0.00000000e+00,
 -1.58920825e+06,  0.00000000e+00,  0.00000000e+00,
  5.24938343e+06,  0.00000000e+00,  0.00000000e+00,
  3.83775006e+06,  0.00000000e+00,  0.00000000e+00,
  1.51006944e+06,  0.00000000e+00,  0.00000000e+00,
  4.92936528e+06,  0.00000000e+00,  0.00000000e+00,
  3.79038495e+06,  0.00000000e+00,  0.00000000e+00,
  4.88200017e+06,  0.00000000e+00,  0.00000000e+00,
  7.20968078e+06,  0.00000000e+00,  0.00000000e+00,
  0.00000000e+00, -1.49447802e+06,  0.00000000e+00,
  0.00000000e+00,  3.06615298e+06,  0.00000000e+00,
  0.00000000e+00, -1.58920825e+06,  0.00000000e+00,
  0.00000000e+00,  5.24938343e+06,  0.00000000e+00,
  0.00000000e+00,  3.83775006e+06,  0.00000000e+00,
  0.00000000e+00,  1.51006944e+06,  0.00000000e+00,
  0.00000000e+00,  4.92936528e+06,  0.00000000e+00,
  0.00000000e+00,  3.79038495e+06,  0.00000000e+00,
  0.00000000e+00,  4.88200017e+06,  0.00000000e+00,
  0.00000000e+00,  7.20968078e+06,  0.00000000e+00,
  0.00000000e+00,  0.00000000e+00, -1.49447802e+06,
  0.00000000e+00,  0.00000000e+00,  3.06615298e+06,
  0.00000000e+00,  0.00000000e+00, -1.58920825e+06,
  0.00000000e+00,  0.00000000e+00,  5.24938343e+06,
  0.00000000e+00,  0.00000000e+00,  3.83775006e+06,
  0.00000000e+00,  0.00000000e+00,  1.51006944e+06,
  0.00000000e+00,  0.00000000e+00,  4.92936528e+06,
  0.00000000e+00,  0.00000000e+00,  3.79038495e+06,
  0.00000000e+00,  0.00000000e+00,  4.88200017e+06,
  0.00000000e+00,  0.00000000e+00,  7.20968078e+06,
};

pylith::feassemble::ElasticityExplicitData3DQuadratic::ElasticityExplicitData3DQuadratic(void)
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

pylith::feassemble::ElasticityExplicitData3DQuadratic::~ElasticityExplicitData3DQuadratic(void)
{}


// End of file
