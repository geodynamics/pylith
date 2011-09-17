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
// This file was generated from python application elasticityexplicitapp.

#include "ElasticityExplicitGravData3DLinear.hh"

const int pylith::feassemble::ElasticityExplicitGravData3DLinear::_spaceDim = 3;

const int pylith::feassemble::ElasticityExplicitGravData3DLinear::_cellDim = 3;

const int pylith::feassemble::ElasticityExplicitGravData3DLinear::_numVertices = 4;

const int pylith::feassemble::ElasticityExplicitGravData3DLinear::_numCells = 1;

const int pylith::feassemble::ElasticityExplicitGravData3DLinear::_numBasis = 4;

const int pylith::feassemble::ElasticityExplicitGravData3DLinear::_numQuadPts = 1;

const char* pylith::feassemble::ElasticityExplicitGravData3DLinear::_matType = "ElasticIsotropic3D";

const char* pylith::feassemble::ElasticityExplicitGravData3DLinear::_matDBFilename = "data/elasticisotropic3d.spatialdb";

const int pylith::feassemble::ElasticityExplicitGravData3DLinear::_matId = 0;

const char* pylith::feassemble::ElasticityExplicitGravData3DLinear::_matLabel = "elastic isotropic 3-D";

const PylithScalar pylith::feassemble::ElasticityExplicitGravData3DLinear::_dt =   1.00000000e-02;

const PylithScalar pylith::feassemble::ElasticityExplicitGravData3DLinear::_gravityVec[] = {
  0.00000000e+00,  0.00000000e+00, -1.00000000e+08,
};

const PylithScalar pylith::feassemble::ElasticityExplicitGravData3DLinear::_vertices[] = {
 -5.00000000e-01, -1.00000000e+00, -5.00000000e-01,
  2.00000000e+00, -5.00000000e-01, -4.00000000e-01,
  1.00000000e+00, -1.00000000e-01, -3.00000000e-01,
 -2.00000000e-01,  5.00000000e-01,  2.00000000e+00,
};

const int pylith::feassemble::ElasticityExplicitGravData3DLinear::_cells[] = {
0,1,2,3,
};

const PylithScalar pylith::feassemble::ElasticityExplicitGravData3DLinear::_verticesRef[] = {
 -1.00000000e+00, -1.00000000e+00, -1.00000000e+00,
  1.00000000e+00, -1.00000000e+00, -1.00000000e+00,
 -1.00000000e+00,  1.00000000e+00, -1.00000000e+00,
 -1.00000000e+00, -1.00000000e+00,  1.00000000e+00,
};

const PylithScalar pylith::feassemble::ElasticityExplicitGravData3DLinear::_quadPts[] = {
 -5.00000000e-01, -5.00000000e-01, -5.00000000e-01,
};

const PylithScalar pylith::feassemble::ElasticityExplicitGravData3DLinear::_quadWts[] = {
  1.33333333e+00,
};

const PylithScalar pylith::feassemble::ElasticityExplicitGravData3DLinear::_basis[] = {
  2.50000000e-01,  2.50000000e-01,  2.50000000e-01,
  2.50000000e-01,};

const PylithScalar pylith::feassemble::ElasticityExplicitGravData3DLinear::_basisDerivRef[] = {
 -5.00000000e-01, -5.00000000e-01, -5.00000000e-01,
  5.00000000e-01,  0.00000000e+00,  0.00000000e+00,
  0.00000000e+00,  5.00000000e-01,  0.00000000e+00,
  0.00000000e+00,  0.00000000e+00,  5.00000000e-01,
};

const PylithScalar pylith::feassemble::ElasticityExplicitGravData3DLinear::_fieldTIncr[] = {
  3.00000000e-01,  2.00000000e-01, -5.00000000e-01,
 -3.00000000e-01, -4.00000000e-01, -6.00000000e-01,
  2.00000000e-01,  6.00000000e-01,  3.00000000e-01,
 -6.00000000e-01, -1.00000000e-01, -3.00000000e-01,
};

const PylithScalar pylith::feassemble::ElasticityExplicitGravData3DLinear::_fieldT[] = {
  8.00000000e-01,  1.00000000e-01, -6.00000000e-01,
 -1.00000000e-01, -2.00000000e-01, -5.00000000e-01,
  1.00000000e-01,  7.00000000e-01,  2.00000000e-01,
 -5.00000000e-01, -0.00000000e+00, -2.00000000e-01,
};

const PylithScalar pylith::feassemble::ElasticityExplicitGravData3DLinear::_fieldTmdt[] = {
  1.00000000e-01,  1.00000000e-01, -3.00000000e-01,
 -2.00000000e-01, -1.00000000e-01, -5.00000000e-01,
  2.00000000e-01,  4.00000000e-01,  1.00000000e-01,
 -4.00000000e-01, -1.00000000e-01, -1.00000000e-01,
};

const PylithScalar pylith::feassemble::ElasticityExplicitGravData3DLinear::_valsResidual[] = {
 -6.53693819e+09,  3.88079833e+10, -3.01595567e+10,
 -4.32000975e+09,  7.13967240e+10, -9.96561003e+09,
  7.21670494e+09, -1.13026998e+11, -6.77007835e+10,
  3.64360549e+09,  2.82229089e+09, -2.66713597e+10,
};

const PylithScalar pylith::feassemble::ElasticityExplicitGravData3DLinear::_valsJacobian[] = {
  8.40625000e+05,  0.00000000e+00,  0.00000000e+00,
  8.40625000e+05,  0.00000000e+00,  0.00000000e+00,
  8.40625000e+05,  0.00000000e+00,  0.00000000e+00,
  8.40625000e+05,  0.00000000e+00,  0.00000000e+00,
  0.00000000e+00,  8.40625000e+05,  0.00000000e+00,
  0.00000000e+00,  8.40625000e+05,  0.00000000e+00,
  0.00000000e+00,  8.40625000e+05,  0.00000000e+00,
  0.00000000e+00,  8.40625000e+05,  0.00000000e+00,
  0.00000000e+00,  0.00000000e+00,  8.40625000e+05,
  0.00000000e+00,  0.00000000e+00,  8.40625000e+05,
  0.00000000e+00,  0.00000000e+00,  8.40625000e+05,
  0.00000000e+00,  0.00000000e+00,  8.40625000e+05,
  8.40625000e+05,  0.00000000e+00,  0.00000000e+00,
  8.40625000e+05,  0.00000000e+00,  0.00000000e+00,
  8.40625000e+05,  0.00000000e+00,  0.00000000e+00,
  8.40625000e+05,  0.00000000e+00,  0.00000000e+00,
  0.00000000e+00,  8.40625000e+05,  0.00000000e+00,
  0.00000000e+00,  8.40625000e+05,  0.00000000e+00,
  0.00000000e+00,  8.40625000e+05,  0.00000000e+00,
  0.00000000e+00,  8.40625000e+05,  0.00000000e+00,
  0.00000000e+00,  0.00000000e+00,  8.40625000e+05,
  0.00000000e+00,  0.00000000e+00,  8.40625000e+05,
  0.00000000e+00,  0.00000000e+00,  8.40625000e+05,
  0.00000000e+00,  0.00000000e+00,  8.40625000e+05,
  8.40625000e+05,  0.00000000e+00,  0.00000000e+00,
  8.40625000e+05,  0.00000000e+00,  0.00000000e+00,
  8.40625000e+05,  0.00000000e+00,  0.00000000e+00,
  8.40625000e+05,  0.00000000e+00,  0.00000000e+00,
  0.00000000e+00,  8.40625000e+05,  0.00000000e+00,
  0.00000000e+00,  8.40625000e+05,  0.00000000e+00,
  0.00000000e+00,  8.40625000e+05,  0.00000000e+00,
  0.00000000e+00,  8.40625000e+05,  0.00000000e+00,
  0.00000000e+00,  0.00000000e+00,  8.40625000e+05,
  0.00000000e+00,  0.00000000e+00,  8.40625000e+05,
  0.00000000e+00,  0.00000000e+00,  8.40625000e+05,
  0.00000000e+00,  0.00000000e+00,  8.40625000e+05,
  8.40625000e+05,  0.00000000e+00,  0.00000000e+00,
  8.40625000e+05,  0.00000000e+00,  0.00000000e+00,
  8.40625000e+05,  0.00000000e+00,  0.00000000e+00,
  8.40625000e+05,  0.00000000e+00,  0.00000000e+00,
  0.00000000e+00,  8.40625000e+05,  0.00000000e+00,
  0.00000000e+00,  8.40625000e+05,  0.00000000e+00,
  0.00000000e+00,  8.40625000e+05,  0.00000000e+00,
  0.00000000e+00,  8.40625000e+05,  0.00000000e+00,
  0.00000000e+00,  0.00000000e+00,  8.40625000e+05,
  0.00000000e+00,  0.00000000e+00,  8.40625000e+05,
  0.00000000e+00,  0.00000000e+00,  8.40625000e+05,
  0.00000000e+00,  0.00000000e+00,  8.40625000e+05,
};

const PylithScalar pylith::feassemble::ElasticityExplicitGravData3DLinear::_valsResidualLumped[] = {
 -6.53643381e+09,  3.88073108e+10, -3.01595567e+10,
 -4.31950537e+09,  7.13977327e+10, -9.96426503e+09,
  7.21485556e+09, -1.13028007e+11, -6.77021285e+10,
  3.64444612e+09,  2.82296339e+09, -2.66713597e+10,
};

const PylithScalar pylith::feassemble::ElasticityExplicitGravData3DLinear::_valsJacobianLumped[] = {
  3.36250000e+06,  3.36250000e+06,  3.36250000e+06,
  3.36250000e+06,  3.36250000e+06,  3.36250000e+06,
  3.36250000e+06,  3.36250000e+06,  3.36250000e+06,
  3.36250000e+06,  3.36250000e+06,  3.36250000e+06,
};

pylith::feassemble::ElasticityExplicitGravData3DLinear::ElasticityExplicitGravData3DLinear(void)
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
  gravityVec = const_cast<PylithScalar*>(_gravityVec);
  vertices = const_cast<PylithScalar*>(_vertices);
  cells = const_cast<int*>(_cells);
  verticesRef = const_cast<PylithScalar*>(_verticesRef);
  quadPts = const_cast<PylithScalar*>(_quadPts);
  quadWts = const_cast<PylithScalar*>(_quadWts);
  basis = const_cast<PylithScalar*>(_basis);
  basisDerivRef = const_cast<PylithScalar*>(_basisDerivRef);
  fieldTIncr = const_cast<PylithScalar*>(_fieldTIncr);
  fieldT = const_cast<PylithScalar*>(_fieldT);
  fieldTmdt = const_cast<PylithScalar*>(_fieldTmdt);
  valsResidual = const_cast<PylithScalar*>(_valsResidual);
  valsJacobian = const_cast<PylithScalar*>(_valsJacobian);
  valsResidualLumped = const_cast<PylithScalar*>(_valsResidualLumped);
  valsJacobianLumped = const_cast<PylithScalar*>(_valsJacobianLumped);
} // constructor

pylith::feassemble::ElasticityExplicitGravData3DLinear::~ElasticityExplicitGravData3DLinear(void)
{}


// End of file
