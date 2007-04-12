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
// This file was generated from python application quadrature1din3dquadratic.

#include "QuadratureData1Din3DQuadratic.hh"

const int pylith::feassemble::QuadratureData1Din3DQuadratic::_numVertices = 3;

const int pylith::feassemble::QuadratureData1Din3DQuadratic::_spaceDim = 3;

const int pylith::feassemble::QuadratureData1Din3DQuadratic::_numCells = 1;

const int pylith::feassemble::QuadratureData1Din3DQuadratic::_cellDim = 1;

const int pylith::feassemble::QuadratureData1Din3DQuadratic::_numBasis = 3;

const int pylith::feassemble::QuadratureData1Din3DQuadratic::_numQuadPts = 2;

const double pylith::feassemble::QuadratureData1Din3DQuadratic::_vertices[] = {
  1.00000000e+00, -1.50000000e+00, -2.00000000e+00,
  3.00000000e-01,  3.00000000e-01,  8.00000000e-01,
 -5.00000000e-01,  2.00000000e+00,  3.00000000e+00,
};

const int pylith::feassemble::QuadratureData1Din3DQuadratic::_cells[] = {
       0,       1,       2,
};

const double pylith::feassemble::QuadratureData1Din3DQuadratic::_quadPtsRef[] = {
 -5.77350269e-01,
  5.77350269e-01,
};

const double pylith::feassemble::QuadratureData1Din3DQuadratic::_quadWts[] = {
  1.00000000e+00,  1.00000000e+00,
};

const double pylith::feassemble::QuadratureData1Din3DQuadratic::_basis[] = {
  4.55341801e-01,
  6.66666667e-01,
 -1.22008468e-01,
 -1.22008468e-01,
  6.66666667e-01,
  4.55341801e-01,
};

const double pylith::feassemble::QuadratureData1Din3DQuadratic::_basisDeriv[] = {
 -1.07735027e+00,
  1.15470054e+00,
 -7.73502692e-02,
  7.73502692e-02,
 -1.15470054e+00,
  1.07735027e+00,
};

const double pylith::feassemble::QuadratureData1Din3DQuadratic::_quadPts[] = {
  7.16346035e-01, -7.27029638e-01, -7.43375673e-01,
 -1.49679369e-01,  1.29369630e+00,  2.14337567e+00,
};

const double pylith::feassemble::QuadratureData1Din3DQuadratic::_jacobian[] = {
 -6.92264973e-01,  1.80773503e+00,  2.84641016e+00,
 -8.07735027e-01,  1.69226497e+00,  2.15358984e+00,
};

const double pylith::feassemble::QuadratureData1Din3DQuadratic::_jacobianDet[] = {
  3.44226488e+00,  2.85554650e+00,
};

const double pylith::feassemble::QuadratureData1Din3DQuadratic::_jacobianInv[] = {
 -1.44453358e+00,
  5.53178417e-01,
  3.51319713e-01,
 -1.23802976e+00,
  5.90924008e-01,
  4.64340973e-01,
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
  quadPtsRef = const_cast<double*>(_quadPtsRef);
  quadWts = const_cast<double*>(_quadWts);
  basis = const_cast<double*>(_basis);
  basisDeriv = const_cast<double*>(_basisDeriv);
  quadPts = const_cast<double*>(_quadPts);
  jacobian = const_cast<double*>(_jacobian);
  jacobianDet = const_cast<double*>(_jacobianDet);
  jacobianInv = const_cast<double*>(_jacobianInv);
} // constructor

pylith::feassemble::QuadratureData1Din3DQuadratic::~QuadratureData1Din3DQuadratic(void)
{}


// End of file
