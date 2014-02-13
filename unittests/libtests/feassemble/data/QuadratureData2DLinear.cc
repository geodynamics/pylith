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
// Copyright (c) 2010-2014 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

// DO NOT EDIT THIS FILE
// This file was generated from python application quadratureapp.

#include "QuadratureData2DLinear.hh"

const int pylith::feassemble::QuadratureData2DLinear::_numVertices = 3;

const int pylith::feassemble::QuadratureData2DLinear::_spaceDim = 2;

const int pylith::feassemble::QuadratureData2DLinear::_numCells = 1;

const int pylith::feassemble::QuadratureData2DLinear::_cellDim = 2;

const int pylith::feassemble::QuadratureData2DLinear::_numBasis = 3;

const int pylith::feassemble::QuadratureData2DLinear::_numQuadPts = 1;

const PylithScalar pylith::feassemble::QuadratureData2DLinear::_vertices[] = {
  2.00000000e-01, -4.00000000e-01,
  3.00000000e-01,  5.00000000e-01,
 -1.00000000e+00, -2.00000000e-01,
};

const int pylith::feassemble::QuadratureData2DLinear::_cells[] = {
       0,       1,       2,
};

const PylithScalar pylith::feassemble::QuadratureData2DLinear::_verticesRef[] = {
 -1.00000000e+00, -1.00000000e+00,
  1.00000000e+00, -1.00000000e+00,
 -1.00000000e+00,  1.00000000e+00,
};

const PylithScalar pylith::feassemble::QuadratureData2DLinear::_quadPtsRef[] = {
 -3.33333333e-01, -3.33333333e-01,
};

const PylithScalar pylith::feassemble::QuadratureData2DLinear::_quadWts[] = {
  2.00000000e+00,
};

const PylithScalar pylith::feassemble::QuadratureData2DLinear::_quadPts[] = {
 -1.66666667e-01, -3.33333333e-02,
};

const PylithScalar pylith::feassemble::QuadratureData2DLinear::_basis[] = {
  3.33333333e-01,  3.33333333e-01,
  3.33333333e-01,};

const PylithScalar pylith::feassemble::QuadratureData2DLinear::_basisDerivRef[] = {
 -5.00000000e-01, -5.00000000e-01,
  5.00000000e-01,  0.00000000e+00,
  0.00000000e+00,  5.00000000e-01,
};

const PylithScalar pylith::feassemble::QuadratureData2DLinear::_basisDeriv[] = {
  6.36363636e-01, -1.18181818e+00,
  1.81818182e-01,  1.09090909e+00,
 -8.18181818e-01,  9.09090909e-02,
};

const PylithScalar pylith::feassemble::QuadratureData2DLinear::_jacobian[] = {
  5.00000000e-02, -6.00000000e-01,
  4.50000000e-01,  1.00000000e-01,
};

const PylithScalar pylith::feassemble::QuadratureData2DLinear::_jacobianDet[] = {
  2.75000000e-01,
};

const PylithScalar pylith::feassemble::QuadratureData2DLinear::_jacobianInv[] = {
  3.63636364e-01,  2.18181818e+00,
 -1.63636364e+00,  1.81818182e-01,
};

pylith::feassemble::QuadratureData2DLinear::QuadratureData2DLinear(void)
{ // constructor
  numVertices = _numVertices;
  spaceDim = _spaceDim;
  numCells = _numCells;
  cellDim = _cellDim;
  numBasis = _numBasis;
  numQuadPts = _numQuadPts;
  vertices = const_cast<PylithScalar*>(_vertices);
  cells = const_cast<int*>(_cells);
  verticesRef = const_cast<PylithScalar*>(_verticesRef);
  quadPtsRef = const_cast<PylithScalar*>(_quadPtsRef);
  quadWts = const_cast<PylithScalar*>(_quadWts);
  quadPts = const_cast<PylithScalar*>(_quadPts);
  basis = const_cast<PylithScalar*>(_basis);
  basisDerivRef = const_cast<PylithScalar*>(_basisDerivRef);
  basisDeriv = const_cast<PylithScalar*>(_basisDeriv);
  jacobian = const_cast<PylithScalar*>(_jacobian);
  jacobianDet = const_cast<PylithScalar*>(_jacobianDet);
  jacobianInv = const_cast<PylithScalar*>(_jacobianInv);
} // constructor

pylith::feassemble::QuadratureData2DLinear::~QuadratureData2DLinear(void)
{}


// End of file
