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

#include "QuadratureData1DQuadratic.hh"

const int pylith::feassemble::QuadratureData1DQuadratic::_numVertices = 3;

const int pylith::feassemble::QuadratureData1DQuadratic::_spaceDim = 1;

const int pylith::feassemble::QuadratureData1DQuadratic::_numCells = 1;

const int pylith::feassemble::QuadratureData1DQuadratic::_cellDim = 1;

const int pylith::feassemble::QuadratureData1DQuadratic::_numBasis = 3;

const int pylith::feassemble::QuadratureData1DQuadratic::_numQuadPts = 2;

const double pylith::feassemble::QuadratureData1DQuadratic::_vertices[] = {
 -2.50000000e-01,
  2.00000000e+00,
  8.75000000e-01,
};

const int pylith::feassemble::QuadratureData1DQuadratic::_cells[] = {
       0,       1,       2,
};

const double pylith::feassemble::QuadratureData1DQuadratic::_verticesRef[] = {
 -1.00000000e+00,
  1.00000000e+00,
  0.00000000e+00,
};

const double pylith::feassemble::QuadratureData1DQuadratic::_quadPtsRef[] = {
 -5.77350269e-01,
  5.77350269e-01,
};

const double pylith::feassemble::QuadratureData1DQuadratic::_quadWts[] = {
  1.00000000e+00,  1.00000000e+00,
};

const double pylith::feassemble::QuadratureData1DQuadratic::_quadPts[] = {
  2.25480947e-01,
  1.52451905e+00,
};

const double pylith::feassemble::QuadratureData1DQuadratic::_basis[] = {
  4.55341801e-01,
 -1.22008468e-01,
  6.66666667e-01,
 -1.22008468e-01,
  4.55341801e-01,
  6.66666667e-01,
};

const double pylith::feassemble::QuadratureData1DQuadratic::_basisDerivRef[] = {
 -1.07735027e+00,
 -7.73502692e-02,
  1.15470054e+00,
  7.73502692e-02,
  1.07735027e+00,
 -1.15470054e+00,
};

const double pylith::feassemble::QuadratureData1DQuadratic::_basisDeriv[] = {
 -9.57644684e-01,
 -6.87557948e-02,
  1.02640048e+00,
  6.87557948e-02,
  9.57644684e-01,
 -1.02640048e+00,
};

const double pylith::feassemble::QuadratureData1DQuadratic::_jacobian[] = {
  1.12500000e+00,
  1.12500000e+00,
};

const double pylith::feassemble::QuadratureData1DQuadratic::_jacobianDet[] = {
  1.12500000e+00,  1.12500000e+00,
};

const double pylith::feassemble::QuadratureData1DQuadratic::_jacobianInv[] = {
  8.88888889e-01,
  8.88888889e-01,
};

pylith::feassemble::QuadratureData1DQuadratic::QuadratureData1DQuadratic(void)
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

pylith::feassemble::QuadratureData1DQuadratic::~QuadratureData1DQuadratic(void)
{}


// End of file
