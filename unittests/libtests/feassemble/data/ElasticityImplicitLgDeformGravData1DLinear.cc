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
// Copyright (c) 2010 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

// DO NOT EDIT THIS FILE
// This file was generated from python application elasticitylgdeformapp.

#include "ElasticityImplicitLgDeformGravData1DLinear.hh"

const int pylith::feassemble::ElasticityImplicitLgDeformGravData1DLinear::_spaceDim = 1;

const int pylith::feassemble::ElasticityImplicitLgDeformGravData1DLinear::_cellDim = 1;

const int pylith::feassemble::ElasticityImplicitLgDeformGravData1DLinear::_numVertices = 2;

const int pylith::feassemble::ElasticityImplicitLgDeformGravData1DLinear::_numCells = 1;

const int pylith::feassemble::ElasticityImplicitLgDeformGravData1DLinear::_numBasis = 2;

const int pylith::feassemble::ElasticityImplicitLgDeformGravData1DLinear::_numQuadPts = 1;

const char* pylith::feassemble::ElasticityImplicitLgDeformGravData1DLinear::_matType = "ElasticStrain1D";

const char* pylith::feassemble::ElasticityImplicitLgDeformGravData1DLinear::_matDBFilename = "data/elasticstrain1d.spatialdb";

const int pylith::feassemble::ElasticityImplicitLgDeformGravData1DLinear::_matId = 0;

const char* pylith::feassemble::ElasticityImplicitLgDeformGravData1DLinear::_matLabel = "elastic strain 1-D";

const double pylith::feassemble::ElasticityImplicitLgDeformGravData1DLinear::_dt =   1.00000000e-02;

const double pylith::feassemble::ElasticityImplicitLgDeformGravData1DLinear::_gravityVec[] = {
 -1.00000000e+08,
};

const double pylith::feassemble::ElasticityImplicitLgDeformGravData1DLinear::_vertices[] = {
 -2.50000000e-01,
  2.00000000e+00,
};

const int pylith::feassemble::ElasticityImplicitLgDeformGravData1DLinear::_cells[] = {
0,1,
};

const double pylith::feassemble::ElasticityImplicitLgDeformGravData1DLinear::_verticesRef[] = {
 -1.00000000e+00,
  1.00000000e+00,
};

const double pylith::feassemble::ElasticityImplicitLgDeformGravData1DLinear::_quadPts[] = {
  0.00000000e+00,
};

const double pylith::feassemble::ElasticityImplicitLgDeformGravData1DLinear::_quadWts[] = {
  2.00000000e+00,
};

const double pylith::feassemble::ElasticityImplicitLgDeformGravData1DLinear::_basis[] = {
  5.00000000e-01,
  5.00000000e-01,
};

const double pylith::feassemble::ElasticityImplicitLgDeformGravData1DLinear::_basisDerivRef[] = {
 -5.00000000e-01,
  5.00000000e-01,
};

const double pylith::feassemble::ElasticityImplicitLgDeformGravData1DLinear::_fieldTIncr[] = {
  1.20000000e+00,
  1.70000000e+00,
};

const double pylith::feassemble::ElasticityImplicitLgDeformGravData1DLinear::_fieldT[] = {
  1.10000000e+00,
  1.50000000e+00,
};

const double pylith::feassemble::ElasticityImplicitLgDeformGravData1DLinear::_fieldTmdt[] = {
  1.00000000e+00,
  1.30000000e+00,
};

const double pylith::feassemble::ElasticityImplicitLgDeformGravData1DLinear::_valsResidual[] = {
 -2.20770000e+11,
 -3.41730000e+11,
};

const double pylith::feassemble::ElasticityImplicitLgDeformGravData1DLinear::_valsJacobian[] = {
  9.76000000e+10,
 -9.76000000e+10,
 -9.76000000e+10,
  9.76000000e+10,
};

pylith::feassemble::ElasticityImplicitLgDeformGravData1DLinear::ElasticityImplicitLgDeformGravData1DLinear(void)
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

pylith::feassemble::ElasticityImplicitLgDeformGravData1DLinear::~ElasticityImplicitLgDeformGravData1DLinear(void)
{}


// End of file
