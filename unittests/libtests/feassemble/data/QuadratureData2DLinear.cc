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
// This file was generated from python application quadrature1dlinear.

#include "QuadratureData2DLinear.hh"

const int pylith::feassemble::QuadratureData2DLinear::_numVertices = 3;

const int pylith::feassemble::QuadratureData2DLinear::_spaceDim = 2;

const int pylith::feassemble::QuadratureData2DLinear::_numCells = 1;

const int pylith::feassemble::QuadratureData2DLinear::_cellDim = 2;

const int pylith::feassemble::QuadratureData2DLinear::_numBasis = 3;

const int pylith::feassemble::QuadratureData2DLinear::_numQuadPts = 1;

const double pylith::feassemble::QuadratureData2DLinear::_vertices[] = {
  2.00000000e-01, -4.00000000e-01,
  3.00000000e-01,  5.00000000e-01,
 -1.00000000e+00, -2.00000000e-01,
};

const int pylith::feassemble::QuadratureData2DLinear::_cells[] = {
       0,       1,       2,
};

const double pylith::feassemble::QuadratureData2DLinear::_quadPtsRef[] = {
  3.33333333e-01,  3.33333333e-01,
};

const double pylith::feassemble::QuadratureData2DLinear::_quadWts[] = {
  5.00000000e-01,
};

const double pylith::feassemble::QuadratureData2DLinear::_quadPts[] = {
 -1.66666667e-01, -3.33333333e-02,
};

const double pylith::feassemble::QuadratureData2DLinear::_basisVert[] = {
  1.00000000e+00,  0.00000000e+00,
  0.00000000e+00,  0.00000000e+00,
  1.00000000e+00,  0.00000000e+00,
  0.00000000e+00,  0.00000000e+00,
  1.00000000e+00,};

const double pylith::feassemble::QuadratureData2DLinear::_basisDerivVert[] = {
 -1.00000000e+00, -1.00000000e+00,
  1.00000000e+00,  0.00000000e+00,
  0.00000000e+00,  1.00000000e+00,
 -1.00000000e+00, -1.00000000e+00,
  1.00000000e+00,  0.00000000e+00,
  0.00000000e+00,  1.00000000e+00,
 -1.00000000e+00, -1.00000000e+00,
  1.00000000e+00,  0.00000000e+00,
  0.00000000e+00,  1.00000000e+00,
};

const double pylith::feassemble::QuadratureData2DLinear::_jacobianVert[] = {
  1.00000000e-01,  9.00000000e-01,
 -1.20000000e+00,  2.00000000e-01,
  1.00000000e-01,  9.00000000e-01,
 -1.20000000e+00,  2.00000000e-01,
  1.00000000e-01,  9.00000000e-01,
 -1.20000000e+00,  2.00000000e-01,
};

const double pylith::feassemble::QuadratureData2DLinear::_jacobianDetVert[] = {
  1.10000000e+00,  1.10000000e+00,  1.10000000e+00,
};

const double pylith::feassemble::QuadratureData2DLinear::_basisQuad[] = {
  3.33333333e-01,  3.33333333e-01,
  3.33333333e-01,};

const double pylith::feassemble::QuadratureData2DLinear::_basisDerivQuad[] = {
 -1.00000000e+00, -1.00000000e+00,
  1.00000000e+00,  0.00000000e+00,
  0.00000000e+00,  1.00000000e+00,
};

const double pylith::feassemble::QuadratureData2DLinear::_jacobianQuad[] = {
  1.00000000e-01,  9.00000000e-01,
 -1.20000000e+00,  2.00000000e-01,
};

const double pylith::feassemble::QuadratureData2DLinear::_jacobianDetQuad[] = {
  1.10000000e+00,
};

const double pylith::feassemble::QuadratureData2DLinear::_jacobianInvQuad[] = {
  1.81818182e-01, -8.18181818e-01,
  1.09090909e+00,  9.09090909e-02,
};

pylith::feassemble::QuadratureData2DLinear::QuadratureData2DLinear(void)
{ // constructor
  numVertices = _numVertices;
  spaceDim = _spaceDim;
  numCells = _numCells;
  cellDim = _cellDim;
  numBasis = _numBasis;
  numQuadPts = _numQuadPts;
  vertices = const_cast<double*>(_vertices);
  cells = const_cast<int*>(_cells);
  quadPtsRef = const_cast<double*>(_quadPtsRef);
  quadWts = const_cast<double*>(_quadWts);
  quadPts = const_cast<double*>(_quadPts);
  basisVert = const_cast<double*>(_basisVert);
  basisDerivVert = const_cast<double*>(_basisDerivVert);
  jacobianVert = const_cast<double*>(_jacobianVert);
  jacobianDetVert = const_cast<double*>(_jacobianDetVert);
  basisQuad = const_cast<double*>(_basisQuad);
  basisDerivQuad = const_cast<double*>(_basisDerivQuad);
  jacobianQuad = const_cast<double*>(_jacobianQuad);
  jacobianDetQuad = const_cast<double*>(_jacobianDetQuad);
  jacobianInvQuad = const_cast<double*>(_jacobianInvQuad);
} // constructor

pylith::feassemble::QuadratureData2DLinear::~QuadratureData2DLinear(void)
{}


// End of file
