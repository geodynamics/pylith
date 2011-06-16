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
// This file was generated from python application quadratureapp.

#include "QuadratureData1Din3DQuadratic.hh"

const int pylith::feassemble::QuadratureData1Din3DQuadratic::_numVertices = 3;

const int pylith::feassemble::QuadratureData1Din3DQuadratic::_spaceDim = 3;

const int pylith::feassemble::QuadratureData1Din3DQuadratic::_numCells = 1;

const int pylith::feassemble::QuadratureData1Din3DQuadratic::_cellDim = 1;

const int pylith::feassemble::QuadratureData1Din3DQuadratic::_numBasis = 3;

const int pylith::feassemble::QuadratureData1Din3DQuadratic::_numQuadPts = 2;

const double pylith::feassemble::QuadratureData1Din3DQuadratic::_vertices[] = {
  1.00000000e+00, -1.50000000e+00, -2.00000000e+00,
 -5.00000000e-01,  2.00000000e+00,  3.00000000e+00,
  2.50000000e-01,  2.50000000e-01,  5.00000000e-01,
};

const int pylith::feassemble::QuadratureData1Din3DQuadratic::_cells[] = {
       0,       1,       2,
};

const double pylith::feassemble::QuadratureData1Din3DQuadratic::_verticesRef[] = {
 -1.00000000e+00,
  1.00000000e+00,
  0.00000000e+00,
};

const double pylith::feassemble::QuadratureData1Din3DQuadratic::_quadPtsRef[] = {
 -5.77350269e-01,
  5.77350269e-01,
};

const double pylith::feassemble::QuadratureData1Din3DQuadratic::_quadWts[] = {
  1.00000000e+00,  1.00000000e+00,
};

const double pylith::feassemble::QuadratureData1Din3DQuadratic::_quadPts[] = {
  6.83012702e-01, -7.60362971e-01, -9.43375673e-01,
 -1.83012702e-01,  1.26036297e+00,  1.94337567e+00,
};

const double pylith::feassemble::QuadratureData1Din3DQuadratic::_basis[] = {
  4.55341801e-01,
 -1.22008468e-01,
  6.66666667e-01,
 -1.22008468e-01,
  4.55341801e-01,
  6.66666667e-01,
};

const double pylith::feassemble::QuadratureData1Din3DQuadratic::_basisDerivRef[] = {
 -1.07735027e+00,
 -7.73502692e-02,
  1.15470054e+00,
  7.73502692e-02,
  1.07735027e+00,
 -1.15470054e+00,
};

const double pylith::feassemble::QuadratureData1Din3DQuadratic::_basisDeriv[] = {
  1.43646703e+00, -6.15628725e-01, -4.30940108e-01,
  1.03133692e-01, -4.42001538e-02, -3.09401077e-02,
 -1.53960072e+00,  6.59828879e-01,  4.61880215e-01,
 -1.03133692e-01,  4.42001538e-02,  3.09401077e-02,
 -1.43646703e+00,  6.15628725e-01,  4.30940108e-01,
  1.53960072e+00, -6.59828879e-01, -4.61880215e-01,
};

const double pylith::feassemble::QuadratureData1Din3DQuadratic::_jacobian[] = {
 -7.50000000e-01,
  1.75000000e+00,
  2.50000000e+00,
 -7.50000000e-01,
  1.75000000e+00,
  2.50000000e+00,
};

const double pylith::feassemble::QuadratureData1Din3DQuadratic::_jacobianDet[] = {
  3.14245127e+00,  3.14245127e+00,
};

const double pylith::feassemble::QuadratureData1Din3DQuadratic::_jacobianInv[] = {
 -1.33333333e+00,  5.71428571e-01,  4.00000000e-01,
 -1.33333333e+00,  5.71428571e-01,  4.00000000e-01,
};

pylith::feassemble::QuadratureData1Din3DQuadratic::QuadratureData1Din3DQuadratic(void)
{ // constructor
  numVertices = _numVertices;
  spaceDim = _spaceDim;
  numCells = _numCells;
  cellDim = _cellDim;
  numBasis = _numBasis;
  numQuadPts = _numQuadPts;
  vertices = const_cast<double*>(_vertices);
  cells = const_cast<int*>(_cells);
  verticesRef = const_cast<double*>(_verticesRef);
  quadPtsRef = const_cast<double*>(_quadPtsRef);
  quadWts = const_cast<double*>(_quadWts);
  quadPts = const_cast<double*>(_quadPts);
  basis = const_cast<double*>(_basis);
  basisDerivRef = const_cast<double*>(_basisDerivRef);
  basisDeriv = const_cast<double*>(_basisDeriv);
  jacobian = const_cast<double*>(_jacobian);
  jacobianDet = const_cast<double*>(_jacobianDet);
  jacobianInv = const_cast<double*>(_jacobianInv);
} // constructor

pylith::feassemble::QuadratureData1Din3DQuadratic::~QuadratureData1Din3DQuadratic(void)
{}


// End of file
