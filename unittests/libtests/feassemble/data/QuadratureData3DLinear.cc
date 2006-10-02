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

#include "QuadratureData3DLinear.hh"

const int pylith::feassemble::QuadratureData3DLinear::_numVertices = 4;

const int pylith::feassemble::QuadratureData3DLinear::_spaceDim = 3;

const int pylith::feassemble::QuadratureData3DLinear::_numCells = 1;

const int pylith::feassemble::QuadratureData3DLinear::_cellDim = 3;

const int pylith::feassemble::QuadratureData3DLinear::_numCorners = 4;

const int pylith::feassemble::QuadratureData3DLinear::_numQuadPts = 1;

const double pylith::feassemble::QuadratureData3DLinear::_vertices[] = {
 -5.00000000e-01, -1.00000000e+00, -5.00000000e-01,
  2.00000000e+00, -5.00000000e-01, -4.00000000e-01,
  1.00000000e+00, -1.00000000e-01, -3.00000000e-01,
 -2.00000000e-01,  5.00000000e-01,  2.00000000e+00,
};

const int pylith::feassemble::QuadratureData3DLinear::_cells[] = {
       0,       1,       2,       3,
};

const double pylith::feassemble::QuadratureData3DLinear::_quadPtsRef[] = {
  3.33333333e-01,  3.33333333e-01,  3.33333333e-01,
};

const double pylith::feassemble::QuadratureData3DLinear::_quadWts[] = {
  1.66666667e-01,
};

const double pylith::feassemble::QuadratureData3DLinear::_basis[] = {
  1.11022302e-16,  3.33333333e-01,  3.33333333e-01,
  3.33333333e-01,};

const double pylith::feassemble::QuadratureData3DLinear::_basisDeriv[] = {
 -1.00000000e+00, -1.00000000e+00, -1.00000000e+00,
  1.00000000e+00,  0.00000000e+00,  0.00000000e+00,
  0.00000000e+00,  1.00000000e+00,  0.00000000e+00,
  0.00000000e+00,  0.00000000e+00,  1.00000000e+00,
};

const double pylith::feassemble::QuadratureData3DLinear::_quadPts[] = {
  9.33333333e-01, -3.33333333e-02,  4.33333333e-01,
};

const double pylith::feassemble::QuadratureData3DLinear::_jacobian[] = {
  2.50000000e+00,  5.00000000e-01,  1.00000000e-01,
  1.50000000e+00,  9.00000000e-01,  2.00000000e-01,
  3.00000000e-01,  1.50000000e+00,  2.50000000e+00,
};

const double pylith::feassemble::QuadratureData3DLinear::_jacobianDet[] = {
  3.22800000e+00,
};

const double pylith::feassemble::QuadratureData3DLinear::_jacobianInv[] = {
  6.04089219e-01, -3.40768278e-01,  3.09789343e-03,
 -1.14312268e+00,  1.92688971e+00, -1.08426270e-01,
  6.13382900e-01, -1.11524164e+00,  4.64684015e-01,
};

pylith::feassemble::QuadratureData3DLinear::QuadratureData3DLinear(void)
{ // constructor
  numVertices = _numVertices;
  spaceDim = _spaceDim;
  numCells = _numCells;
  cellDim = _cellDim;
  numCorners = _numCorners;
  numQuadPts = _numQuadPts;
  vertices = const_cast<double*>(_vertices);
  cells = const_cast<int*>(_cells);
  quadPtsRef = const_cast<double*>(_quadPtsRef);
  quadWts = const_cast<double*>(_quadWts);
  basis = const_cast<double*>(_basis);
  basisDeriv = const_cast<double*>(_basisDeriv);
  quadPts = const_cast<double*>(_quadPts);
  jacobian = const_cast<double*>(_jacobian);
  jacobianDet = const_cast<double*>(_jacobianDet);
  jacobianInv = const_cast<double*>(_jacobianInv);
} // constructor

pylith::feassemble::QuadratureData3DLinear::~QuadratureData3DLinear(void)
{}


// End of file
