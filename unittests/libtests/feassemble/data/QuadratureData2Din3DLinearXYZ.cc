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

#include "QuadratureData2Din3DLinearXYZ.hh"

const int pylith::feassemble::QuadratureData2Din3DLinearXYZ::_numVertices = 3;

const int pylith::feassemble::QuadratureData2Din3DLinearXYZ::_spaceDim = 3;

const int pylith::feassemble::QuadratureData2Din3DLinearXYZ::_numCells = 1;

const int pylith::feassemble::QuadratureData2Din3DLinearXYZ::_cellDim = 2;

const int pylith::feassemble::QuadratureData2Din3DLinearXYZ::_numBasis = 3;

const int pylith::feassemble::QuadratureData2Din3DLinearXYZ::_numQuadPts = 1;

const double pylith::feassemble::QuadratureData2Din3DLinearXYZ::_vertices[] = {
  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,
 -1.00000000e+00,  0.00000000e+00,  0.00000000e+00,
  0.00000000e+00,  0.00000000e+00,  1.00000000e+00,
};

const int pylith::feassemble::QuadratureData2Din3DLinearXYZ::_cells[] = {
       0,       1,       2,
};

const double pylith::feassemble::QuadratureData2Din3DLinearXYZ::_verticesRef[] = {
 -1.00000000e+00, -1.00000000e+00,
  1.00000000e+00, -1.00000000e+00,
 -1.00000000e+00,  1.00000000e+00,
};

const double pylith::feassemble::QuadratureData2Din3DLinearXYZ::_quadPtsRef[] = {
  3.33333333e-01,  3.33333333e-01,
};

const double pylith::feassemble::QuadratureData2Din3DLinearXYZ::_quadWts[] = {
  5.00000000e-01,
};

const double pylith::feassemble::QuadratureData2Din3DLinearXYZ::_quadPts[] = {
 -3.33333333e-01,  0.00000000e+00,  3.33333333e-01,
};

const double pylith::feassemble::QuadratureData2Din3DLinearXYZ::_basis[] = {
  3.33333333e-01,  3.33333333e-01,
  3.33333333e-01,};

const double pylith::feassemble::QuadratureData2Din3DLinearXYZ::_basisDerivRef[] = {
 -1.00000000e+00, -1.00000000e+00,
  1.00000000e+00,  0.00000000e+00,
  0.00000000e+00,  1.00000000e+00,
};

const double pylith::feassemble::QuadratureData2Din3DLinearXYZ::_basisDeriv[] = {
  1.00000000e+00, -0.00000000e+00, -1.00000000e+00,
 -1.00000000e+00,  0.00000000e+00,  0.00000000e+00,
  0.00000000e+00,  0.00000000e+00,  1.00000000e+00,
};

const double pylith::feassemble::QuadratureData2Din3DLinearXYZ::_jacobian[] = {
 -1.00000000e+00,  0.00000000e+00,
  0.00000000e+00,  0.00000000e+00,
  0.00000000e+00,  1.00000000e+00,
};

const double pylith::feassemble::QuadratureData2Din3DLinearXYZ::_jacobianDet[] = {
  1.00000000e+00,
};

const double pylith::feassemble::QuadratureData2Din3DLinearXYZ::_jacobianInv[] = {
 -1.00000000e+00,  0.00000000e+00, -0.00000000e+00,
  0.00000000e+00,  0.00000000e+00,  1.00000000e+00,
};

pylith::feassemble::QuadratureData2Din3DLinearXYZ::QuadratureData2Din3DLinearXYZ(void)
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

pylith::feassemble::QuadratureData2Din3DLinearXYZ::~QuadratureData2Din3DLinearXYZ(void)
{}


// End of file
