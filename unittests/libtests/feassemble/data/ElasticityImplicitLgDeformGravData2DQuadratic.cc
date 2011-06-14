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
// This file was generated from python application elasticitylgdeformapp.

#include "ElasticityImplicitLgDeformGravData2DQuadratic.hh"

const int pylith::feassemble::ElasticityImplicitLgDeformGravData2DQuadratic::_spaceDim = 2;

const int pylith::feassemble::ElasticityImplicitLgDeformGravData2DQuadratic::_cellDim = 2;

const int pylith::feassemble::ElasticityImplicitLgDeformGravData2DQuadratic::_numVertices = 6;

const int pylith::feassemble::ElasticityImplicitLgDeformGravData2DQuadratic::_numCells = 1;

const int pylith::feassemble::ElasticityImplicitLgDeformGravData2DQuadratic::_numBasis = 6;

const int pylith::feassemble::ElasticityImplicitLgDeformGravData2DQuadratic::_numQuadPts = 6;

const char* pylith::feassemble::ElasticityImplicitLgDeformGravData2DQuadratic::_matType = "ElasticPlaneStrain";

const char* pylith::feassemble::ElasticityImplicitLgDeformGravData2DQuadratic::_matDBFilename = "data/elasticplanestrain.spatialdb";

const int pylith::feassemble::ElasticityImplicitLgDeformGravData2DQuadratic::_matId = 0;

const char* pylith::feassemble::ElasticityImplicitLgDeformGravData2DQuadratic::_matLabel = "elastic strain 2-D";

const double pylith::feassemble::ElasticityImplicitLgDeformGravData2DQuadratic::_dt =   1.00000000e-02;

const double pylith::feassemble::ElasticityImplicitLgDeformGravData2DQuadratic::_gravityVec[] = {
  0.00000000e+00, -1.00000000e+08,
};

const double pylith::feassemble::ElasticityImplicitLgDeformGravData2DQuadratic::_vertices[] = {
 -1.00000000e+00, -1.00000000e+00,
  1.00000000e+00,  2.00000000e-01,
 -1.50000000e+00,  5.00000000e-01,
 -2.50000000e-01,  3.50000000e-01,
 -1.25000000e+00, -2.50000000e-01,
  0.00000000e+00, -4.00000000e-01,
};

const int pylith::feassemble::ElasticityImplicitLgDeformGravData2DQuadratic::_cells[] = {
0,1,2,3,4,5,
};

const double pylith::feassemble::ElasticityImplicitLgDeformGravData2DQuadratic::_verticesRef[] = {
 -1.00000000e+00, -1.00000000e+00,
  1.00000000e+00, -1.00000000e+00,
 -1.00000000e+00,  1.00000000e+00,
  0.00000000e+00,  0.00000000e+00,
 -1.00000000e+00,  0.00000000e+00,
  0.00000000e+00, -1.00000000e+00,
};

const double pylith::feassemble::ElasticityImplicitLgDeformGravData2DQuadratic::_quadPts[] = {
 -7.50000000e-01, -7.50000000e-01,
  7.50000000e-01, -7.50000000e-01,
 -7.50000000e-01,  7.50000000e-01,
  0.00000000e+00, -7.50000000e-01,
 -7.50000000e-01,  0.00000000e+00,
  2.50000000e-01,  2.50000000e-01,
};

const double pylith::feassemble::ElasticityImplicitLgDeformGravData2DQuadratic::_quadWts[] = {
  3.33333333e-01,  3.33333333e-01,  3.33333333e-01,  3.33333333e-01,  3.33333333e-01,  3.33333333e-01,
};

const double pylith::feassemble::ElasticityImplicitLgDeformGravData2DQuadratic::_basis[] = {
  3.75000000e-01, -9.37500000e-02,
 -9.37500000e-02,  6.25000000e-02,
  3.75000000e-01,  3.75000000e-01,
  0.00000000e+00,  6.56250000e-01,
 -9.37500000e-02,  4.37500000e-01,
 -0.00000000e+00, -0.00000000e+00,
  0.00000000e+00, -9.37500000e-02,
  6.56250000e-01,  4.37500000e-01,
 -0.00000000e+00, -0.00000000e+00,
 -9.37500000e-02,  0.00000000e+00,
 -9.37500000e-02,  2.50000000e-01,
  1.87500000e-01,  7.50000000e-01,
 -9.37500000e-02, -9.37500000e-02,
  0.00000000e+00,  2.50000000e-01,
  7.50000000e-01,  1.87500000e-01,
  3.75000000e-01,  1.56250000e-01,
  1.56250000e-01,  1.56250000e+00,
 -6.25000000e-01, -6.25000000e-01,
};

const double pylith::feassemble::ElasticityImplicitLgDeformGravData2DQuadratic::_basisDerivRef[] = {
 -1.00000000e+00, -1.00000000e+00,
 -2.50000000e-01,  0.00000000e+00,
  0.00000000e+00, -2.50000000e-01,
  2.50000000e-01,  2.50000000e-01,
 -2.50000000e-01,  1.25000000e+00,
  1.25000000e+00, -2.50000000e-01,
  5.00000000e-01,  5.00000000e-01,
  1.25000000e+00,  0.00000000e+00,
  0.00000000e+00, -2.50000000e-01,
  2.50000000e-01,  1.75000000e+00,
 -2.50000000e-01, -2.50000000e-01,
 -1.75000000e+00, -1.75000000e+00,
  5.00000000e-01,  5.00000000e-01,
 -2.50000000e-01,  0.00000000e+00,
  0.00000000e+00,  1.25000000e+00,
  1.75000000e+00,  2.50000000e-01,
 -1.75000000e+00, -1.75000000e+00,
 -2.50000000e-01, -2.50000000e-01,
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
  1.00000000e+00,  1.00000000e+00,
  7.50000000e-01,  0.00000000e+00,
  0.00000000e+00,  7.50000000e-01,
  1.25000000e+00,  1.25000000e+00,
 -1.25000000e+00, -1.75000000e+00,
 -1.75000000e+00, -1.25000000e+00,
};

const double pylith::feassemble::ElasticityImplicitLgDeformGravData2DQuadratic::_fieldTIncr[] = {
 -4.00000000e-01, -6.00000000e-01,
  7.00000000e-01,  8.00000000e-01,
  0.00000000e+00,  2.00000000e-01,
 -5.00000000e-01, -4.00000000e-01,
  3.00000000e-01,  9.00000000e-01,
 -3.00000000e-01, -9.00000000e-01,
};

const double pylith::feassemble::ElasticityImplicitLgDeformGravData2DQuadratic::_fieldT[] = {
 -3.00000000e-01, -4.00000000e-01,
  5.00000000e-01,  6.00000000e-01,
  0.00000000e+00,  1.00000000e-01,
 -2.00000000e-01, -3.00000000e-01,
  2.00000000e-01,  3.00000000e-01,
 -1.00000000e-01, -2.00000000e-01,
};

const double pylith::feassemble::ElasticityImplicitLgDeformGravData2DQuadratic::_fieldTmdt[] = {
 -2.00000000e-01, -3.00000000e-01,
  3.00000000e-01,  4.00000000e-01,
  0.00000000e+00, -1.00000000e-01,
 -3.00000000e-01, -2.00000000e-01,
  1.00000000e-01,  4.00000000e-01,
 -2.00000000e-01, -6.00000000e-01,
};

const double pylith::feassemble::ElasticityImplicitLgDeformGravData2DQuadratic::_valsResidual[] = {
  1.03829775e+12,  1.72143634e+12,
 -7.52244241e+11, -8.67453311e+11,
  8.03865887e+11,  1.02236991e+12,
  1.48071362e+12,  5.42069459e+11,
 -2.20475628e+12, -3.86234073e+12,
 -3.65876743e+11,  9.93918331e+11,
};

const double pylith::feassemble::ElasticityImplicitLgDeformGravData2DQuadratic::_valsJacobian[] = {
  1.18059576e+12,  5.56624427e+11,
  2.17020537e+10,  2.82688919e+10,
  7.03170548e+11,  3.46233237e+11,
  9.30815156e+11,  3.73904122e+11,
 -1.97360381e+12, -9.41768992e+11,
 -8.62679703e+11, -3.63261686e+11,
  5.56624427e+11,  1.72905310e+12,
  2.55818801e+10,  7.57022146e+10,
  3.45077232e+11,  7.84996129e+11,
  3.77747140e+11,  6.31931179e+11,
 -9.44649933e+11, -2.68111294e+12,
 -3.60380746e+11, -5.40569680e+11,
  2.17020537e+10,  2.55818801e+10,
  5.70555199e+11,  2.59713039e+11,
 -2.96119462e+10, -2.59165131e+10,
 -4.24956826e+11, -1.36798725e+11,
  1.42337474e+11,  6.73297341e+10,
 -2.80025954e+11, -1.89909415e+11,
  2.82688919e+10,  7.57022146e+10,
  2.59713039e+11,  5.77356324e+11,
 -2.66709077e+10, -5.44779112e+10,
 -1.27273335e+11, -2.04209092e+11,
  6.53971169e+10,  2.26338676e+11,
 -1.99434806e+11, -6.20710212e+11,
  7.03170548e+11,  3.45077232e+11,
 -2.96119462e+10, -2.66709077e+10,
  9.20384196e+11,  4.11776807e+11,
  4.71950241e+11,  3.08031050e+11,
 -1.68844719e+12, -8.05938575e+11,
 -3.77445853e+11, -2.32275605e+11,
  3.46233237e+11,  7.84996129e+11,
 -2.59165131e+10, -5.44779112e+10,
  4.11776807e+11,  9.45887355e+11,
  3.13343794e+11,  3.21137311e+11,
 -8.11251320e+11, -1.83366211e+12,
 -2.34186006e+11, -1.63880774e+11,
  9.30815156e+11,  3.77747140e+11,
 -4.24956826e+11, -1.27273335e+11,
  4.71950241e+11,  3.13343794e+11,
  2.66933376e+12,  8.66541719e+11,
 -1.81944099e+12, -9.28979272e+11,
 -1.82770134e+12, -5.01380046e+11,
  3.73904122e+11,  6.31931179e+11,
 -1.36798725e+11, -2.04209092e+11,
  3.08031050e+11,  3.21137311e+11,
  8.66541719e+11,  2.75153025e+12,
 -9.20742130e+11, -1.64708940e+12,
 -4.90936036e+11, -1.85330024e+12,
 -1.97360381e+12, -9.44649933e+11,
  1.42337474e+11,  6.53971169e+10,
 -1.68844719e+12, -8.11251320e+11,
 -1.81944099e+12, -9.20742130e+11,
  4.26074832e+12,  2.04726998e+12,
  1.07840620e+12,  5.63976291e+11,
 -9.41768992e+11, -2.68111294e+12,
  6.73297341e+10,  2.26338676e+11,
 -8.05938575e+11, -1.83366211e+12,
 -9.28979272e+11, -1.64708940e+12,
  2.04726998e+12,  5.76761534e+12,
  5.62087130e+11,  1.67910439e+11,
 -8.62679703e+11, -3.60380746e+11,
 -2.80025954e+11, -1.99434806e+11,
 -3.77445853e+11, -2.34186006e+11,
 -1.82770134e+12, -4.90936036e+11,
  1.07840620e+12,  5.62087130e+11,
  2.26944665e+12,  7.22850462e+11,
 -3.63261686e+11, -5.40569680e+11,
 -1.89909415e+11, -6.20710212e+11,
 -2.32275605e+11, -1.63880774e+11,
 -5.01380046e+11, -1.85330024e+12,
  5.63976291e+11,  1.67910439e+11,
  7.22850462e+11,  3.01055047e+12,
};

pylith::feassemble::ElasticityImplicitLgDeformGravData2DQuadratic::ElasticityImplicitLgDeformGravData2DQuadratic(void)
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
} // constructor

pylith::feassemble::ElasticityImplicitLgDeformGravData2DQuadratic::~ElasticityImplicitLgDeformGravData2DQuadratic(void)
{}


// End of file
