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

const double pylith::feassemble::ElasticityExplicitData3DQuadratic::_basisDeriv[] = {
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
  1.37303139e+12,  3.59478514e+11,  2.01772066e+12,
 -1.31068344e+11, -1.24490591e+11,  7.74704496e+09,
 -3.68695119e+11, -2.72245362e+12, -1.48103804e+11,
  5.48007226e+11,  3.15715495e+11,  2.06777127e+12,
 -1.21553554e+12, -2.00599538e+12, -3.41366551e+12,
 -1.41823397e+11,  2.78298871e+12, -2.96590241e+11,
 -6.77814032e+12, -6.14553903e+12, -4.48612861e+12,
  4.50987749e+11,  1.95705297e+12,  3.90071124e+11,
  9.75436695e+11,  4.73813145e+12,  2.71900074e+11,
  5.28776613e+12,  8.45145678e+11,  3.58929700e+12,
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
  basisDeriv = const_cast<double*>(_basisDeriv);
  fieldTpdt = const_cast<double*>(_fieldTpdt);
  fieldT = const_cast<double*>(_fieldT);
  fieldTmdt = const_cast<double*>(_fieldTmdt);
  valsResidual = const_cast<double*>(_valsResidual);
  valsJacobian = const_cast<double*>(_valsJacobian);
} // constructor

pylith::feassemble::ElasticityExplicitData3DQuadratic::~ElasticityExplicitData3DQuadratic(void)
{}


// End of file
