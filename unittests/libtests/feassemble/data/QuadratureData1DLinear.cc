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
// This file was generated from python application quadratureapp.

#include "QuadratureData1DLinear.hh"

const int pylith::feassemble::QuadratureData1DLinear::_numVertices = 2;

const int pylith::feassemble::QuadratureData1DLinear::_spaceDim = 1;

const int pylith::feassemble::QuadratureData1DLinear::_numCells = 1;

const int pylith::feassemble::QuadratureData1DLinear::_cellDim = 1;

const int pylith::feassemble::QuadratureData1DLinear::_numBasis = 2;

const int pylith::feassemble::QuadratureData1DLinear::_numQuadPts = 1;

const double pylith::feassemble::QuadratureData1DLinear::_vertices[] = {
 -2.50000000e-01,
  2.00000000e+00,
};

const int pylith::feassemble::QuadratureData1DLinear::_cells[] = {
       0,       1,
};

const double pylith::feassemble::QuadratureData1DLinear::_verticesRef[] = {
 -1.00000000e+00,
  1.00000000e+00,
};

const double pylith::feassemble::QuadratureData1DLinear::_quadPtsRef[] = {
  0.00000000e+00,
};

const double pylith::feassemble::QuadratureData1DLinear::_quadWts[] = {
  2.00000000e+00,
};

const double pylith::feassemble::QuadratureData1DLinear::_quadPts[] = {
  8.75000000e-01,
};

const double pylith::feassemble::QuadratureData1DLinear::_basis[] = {
  5.00000000e-01,
  5.00000000e-01,
};

const double pylith::feassemble::QuadratureData1DLinear::_basisDerivRef[] = {
 -5.00000000e-01,
  5.00000000e-01,
};

const double pylith::feassemble::QuadratureData1DLinear::_basisDeriv[] = {
 -4.44444444e-01,
  4.44444444e-01,
};

const double pylith::feassemble::QuadratureData1DLinear::_jacobian[] = {
  1.12500000e+00,
};

const double pylith::feassemble::QuadratureData1DLinear::_jacobianDet[] = {
  1.12500000e+00,
};

const double pylith::feassemble::QuadratureData1DLinear::_jacobianInv[] = {
  8.88888889e-01,
};

pylith::feassemble::QuadratureData1DLinear::QuadratureData1DLinear(void)
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

pylith::feassemble::QuadratureData1DLinear::~QuadratureData1DLinear(void)
{}


// End of file
