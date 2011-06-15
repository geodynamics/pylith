// -*- C++ -*-
//
// ======================================================================
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2011 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

// DO NOT EDIT THIS FILE
// This file was generated from python application elasticitylgdeformexplicitapp.

#include "ElasticityExplicitLgDeformGravData3DQuadratic.hh"

const int pylith::feassemble::ElasticityExplicitLgDeformGravData3DQuadratic::_spaceDim = 3;

const int pylith::feassemble::ElasticityExplicitLgDeformGravData3DQuadratic::_cellDim = 3;

const int pylith::feassemble::ElasticityExplicitLgDeformGravData3DQuadratic::_numVertices = 10;

const int pylith::feassemble::ElasticityExplicitLgDeformGravData3DQuadratic::_numCells = 1;

const int pylith::feassemble::ElasticityExplicitLgDeformGravData3DQuadratic::_numBasis = 10;

const int pylith::feassemble::ElasticityExplicitLgDeformGravData3DQuadratic::_numQuadPts = 4;

const char* pylith::feassemble::ElasticityExplicitLgDeformGravData3DQuadratic::_matType = "ElasticIsotropic3D";

const char* pylith::feassemble::ElasticityExplicitLgDeformGravData3DQuadratic::_matDBFilename = "data/elasticisotropic3d.spatialdb";

const int pylith::feassemble::ElasticityExplicitLgDeformGravData3DQuadratic::_matId = 0;

const char* pylith::feassemble::ElasticityExplicitLgDeformGravData3DQuadratic::_matLabel = "elastic isotropic 3-D";

const double pylith::feassemble::ElasticityExplicitLgDeformGravData3DQuadratic::_dt =   1.00000000e-02;

const double pylith::feassemble::ElasticityExplicitLgDeformGravData3DQuadratic::_gravityVec[] = {
  0.00000000e+00,  0.00000000e+00, -1.00000000e+08,
};

const double pylith::feassemble::ElasticityExplicitLgDeformGravData3DQuadratic::_vertices[] = {
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

const int pylith::feassemble::ElasticityExplicitLgDeformGravData3DQuadratic::_cells[] = {
0,1,2,3,4,5,6,7,8,9,
};

const double pylith::feassemble::ElasticityExplicitLgDeformGravData3DQuadratic::_verticesRef[] = {
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

const double pylith::feassemble::ElasticityExplicitLgDeformGravData3DQuadratic::_quadPts[] = {
 -8.00000000e-01, -8.00000000e-01, -8.00000000e-01,
  5.00000000e-01, -8.00000000e-01, -8.00000000e-01,
 -8.00000000e-01,  5.00000000e-01, -8.00000000e-01,
 -8.00000000e-01, -8.00000000e-01,  5.00000000e-01,
};

const double pylith::feassemble::ElasticityExplicitLgDeformGravData3DQuadratic::_quadWts[] = {
  3.33333333e-01,  3.33333333e-01,  3.33333333e-01,  3.33333333e-01,
};

const double pylith::feassemble::ElasticityExplicitLgDeformGravData3DQuadratic::_basis[] = {
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

const double pylith::feassemble::ElasticityExplicitLgDeformGravData3DQuadratic::_basisDerivRef[] = {
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

const double pylith::feassemble::ElasticityExplicitLgDeformGravData3DQuadratic::_fieldTIncr[] = {
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

const double pylith::feassemble::ElasticityExplicitLgDeformGravData3DQuadratic::_fieldT[] = {
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

const double pylith::feassemble::ElasticityExplicitLgDeformGravData3DQuadratic::_fieldTmdt[] = {
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

const double pylith::feassemble::ElasticityExplicitLgDeformGravData3DQuadratic::_valsResidual[] = {
  3.68904302e+10,  2.89154635e+10, -2.66246678e+10,
  1.31658958e+11, -3.38188176e+10,  2.73167999e+10,
 -3.92909711e+11, -2.10630391e+11, -2.25786587e+11,
 -1.14269176e+11,  1.30717232e+11, -1.31985028e+11,
 -2.64258575e+11, -2.85535726e+10, -1.61155147e+11,
  4.59398355e+11,  2.32376498e+11,  1.91367337e+11,
 -2.46994434e+11, -7.38524287e+10, -1.31572358e+11,
  5.25219225e+10, -8.35281650e+10, -9.84143036e+10,
 -1.20392699e+11,  3.07708241e+11, -3.84670356e+11,
  4.58368848e+11, -2.69312012e+11,  2.30080394e+11,
};

const double pylith::feassemble::ElasticityExplicitLgDeformGravData3DQuadratic::_valsJacobian[] = {
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

const double pylith::feassemble::ElasticityExplicitLgDeformGravData3DQuadratic::_valsResidualLumped[] = {
  3.68908573e+10,  2.89160572e+10, -2.66286702e+10,
  1.31658131e+11, -3.38201739e+10,  2.73202854e+10,
 -3.92910621e+11, -2.10630330e+11, -2.25786601e+11,
 -1.14267995e+11,  1.30716630e+11, -1.31984968e+11,
 -2.64249877e+11, -2.85444292e+10, -1.61162315e+11,
  4.59400261e+11,  2.32380003e+11,  1.91372671e+11,
 -2.47002605e+11, -7.38605515e+10, -1.31570754e+11,
  5.25288799e+10, -8.35277685e+10, -9.84076887e+10,
 -1.20394352e+11,  3.07717391e+11, -3.84684352e+11,
  4.58361238e+11, -2.69324783e+11,  2.30088475e+11,
};

const double pylith::feassemble::ElasticityExplicitLgDeformGravData3DQuadratic::_valsJacobianLumped[] = {
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

pylith::feassemble::ElasticityExplicitLgDeformGravData3DQuadratic::ElasticityExplicitLgDeformGravData3DQuadratic(void)
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

pylith::feassemble::ElasticityExplicitLgDeformGravData3DQuadratic::~ElasticityExplicitLgDeformGravData3DQuadratic(void)
{}


// End of file
