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
// This file was generated from python application integratorinertia2din3dtwo.

#include "IntegratorDataInertia2Din3DTwo.hh"

const int pylith::feassemble::IntegratorDataInertia2Din3DTwo::_numVertices = 5;

const int pylith::feassemble::IntegratorDataInertia2Din3DTwo::_spaceDim = 3;

const int pylith::feassemble::IntegratorDataInertia2Din3DTwo::_numCells = 2;

const int pylith::feassemble::IntegratorDataInertia2Din3DTwo::_cellDim = 2;

const int pylith::feassemble::IntegratorDataInertia2Din3DTwo::_numCorners = 3;

const int pylith::feassemble::IntegratorDataInertia2Din3DTwo::_numQuadPts = 1;

const int pylith::feassemble::IntegratorDataInertia2Din3DTwo::_fiberDim = 3;

const double pylith::feassemble::IntegratorDataInertia2Din3DTwo::_vertices[] = {
  5.00000000e-01, -2.00000000e+00, -5.00000000e-01,
  3.00000000e+00,  5.00000000e-01,  0.00000000e+00,
 -1.00000000e+00,  2.00000000e+00,  4.00000000e+00,
  7.00000000e+00, -1.00000000e+00, -4.00000000e+00,
  5.50000000e+00,  3.00000000e+00,  5.00000000e-01,
};

const int pylith::feassemble::IntegratorDataInertia2Din3DTwo::_cells[] = {
       0,       1,       2,       3,       4,
       1,};

const double pylith::feassemble::IntegratorDataInertia2Din3DTwo::_quadPts[] = {
  3.33333333e-01,  3.33333333e-01,
};

const double pylith::feassemble::IntegratorDataInertia2Din3DTwo::_quadWts[] = {
  5.00000000e-01,
};

const double pylith::feassemble::IntegratorDataInertia2Din3DTwo::_basis[] = {
  3.33333333e-01,  3.33333333e-01,
  3.33333333e-01,};

const double pylith::feassemble::IntegratorDataInertia2Din3DTwo::_basisDeriv[] = {
 -1.00000000e+00, -1.00000000e+00,
  1.00000000e+00,  0.00000000e+00,
  0.00000000e+00,  1.00000000e+00,
};

const double pylith::feassemble::IntegratorDataInertia2Din3DTwo::_fieldIn[] = {
  1.20000000e+00,  1.00000000e-01, -3.00000000e-01,
  2.00000000e-01, -8.00000000e-01,  1.20000000e+00,
 -9.00000000e-01, -7.00000000e-01,  5.00000000e-01,
  1.30000000e+00, -2.00000000e-01,  1.70000000e+00,
  1.10000000e+00,  1.40000000e+00,  9.00000000e-01,
};

const double pylith::feassemble::IntegratorDataInertia2Din3DTwo::_valsAction[] = {
  5.68342430e-01, -1.59135880e+00,  1.59135880e+00,
  3.52372306e+00, -1.13668486e+00,  5.91076127e+00,
  5.68342430e-01, -1.59135880e+00,  1.59135880e+00,
  2.95538063e+00,  4.54673944e-01,  4.31940246e+00,
  2.95538063e+00,  4.54673944e-01,  4.31940246e+00,
};

const double pylith::feassemble::IntegratorDataInertia2Din3DTwo::_valsMatrix[] = {
  1.13668486e+00,  0.00000000e+00,  0.00000000e+00,
  1.13668486e+00,  0.00000000e+00,  0.00000000e+00,
  1.13668486e+00,  0.00000000e+00,  0.00000000e+00,
  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,
  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,
  0.00000000e+00,  1.13668486e+00,  0.00000000e+00,
  0.00000000e+00,  1.13668486e+00,  0.00000000e+00,
  0.00000000e+00,  1.13668486e+00,  0.00000000e+00,
  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,
  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,
  0.00000000e+00,  0.00000000e+00,  1.13668486e+00,
  0.00000000e+00,  0.00000000e+00,  1.13668486e+00,
  0.00000000e+00,  0.00000000e+00,  1.13668486e+00,
  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,
  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,
  1.13668486e+00,  0.00000000e+00,  0.00000000e+00,
  2.27336972e+00,  0.00000000e+00,  0.00000000e+00,
  1.13668486e+00,  0.00000000e+00,  0.00000000e+00,
  1.13668486e+00,  0.00000000e+00,  0.00000000e+00,
  1.13668486e+00,  0.00000000e+00,  0.00000000e+00,
  0.00000000e+00,  1.13668486e+00,  0.00000000e+00,
  0.00000000e+00,  2.27336972e+00,  0.00000000e+00,
  0.00000000e+00,  1.13668486e+00,  0.00000000e+00,
  0.00000000e+00,  1.13668486e+00,  0.00000000e+00,
  0.00000000e+00,  1.13668486e+00,  0.00000000e+00,
  0.00000000e+00,  0.00000000e+00,  1.13668486e+00,
  0.00000000e+00,  0.00000000e+00,  2.27336972e+00,
  0.00000000e+00,  0.00000000e+00,  1.13668486e+00,
  0.00000000e+00,  0.00000000e+00,  1.13668486e+00,
  0.00000000e+00,  0.00000000e+00,  1.13668486e+00,
  1.13668486e+00,  0.00000000e+00,  0.00000000e+00,
  1.13668486e+00,  0.00000000e+00,  0.00000000e+00,
  1.13668486e+00,  0.00000000e+00,  0.00000000e+00,
  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,
  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,
  0.00000000e+00,  1.13668486e+00,  0.00000000e+00,
  0.00000000e+00,  1.13668486e+00,  0.00000000e+00,
  0.00000000e+00,  1.13668486e+00,  0.00000000e+00,
  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,
  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,
  0.00000000e+00,  0.00000000e+00,  1.13668486e+00,
  0.00000000e+00,  0.00000000e+00,  1.13668486e+00,
  0.00000000e+00,  0.00000000e+00,  1.13668486e+00,
  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,
  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,
  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,
  1.13668486e+00,  0.00000000e+00,  0.00000000e+00,
  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,
  1.13668486e+00,  0.00000000e+00,  0.00000000e+00,
  1.13668486e+00,  0.00000000e+00,  0.00000000e+00,
  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,
  0.00000000e+00,  1.13668486e+00,  0.00000000e+00,
  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,
  0.00000000e+00,  1.13668486e+00,  0.00000000e+00,
  0.00000000e+00,  1.13668486e+00,  0.00000000e+00,
  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,
  0.00000000e+00,  0.00000000e+00,  1.13668486e+00,
  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,
  0.00000000e+00,  0.00000000e+00,  1.13668486e+00,
  0.00000000e+00,  0.00000000e+00,  1.13668486e+00,
  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,
  1.13668486e+00,  0.00000000e+00,  0.00000000e+00,
  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,
  1.13668486e+00,  0.00000000e+00,  0.00000000e+00,
  1.13668486e+00,  0.00000000e+00,  0.00000000e+00,
  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,
  0.00000000e+00,  1.13668486e+00,  0.00000000e+00,
  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,
  0.00000000e+00,  1.13668486e+00,  0.00000000e+00,
  0.00000000e+00,  1.13668486e+00,  0.00000000e+00,
  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,
  0.00000000e+00,  0.00000000e+00,  1.13668486e+00,
  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,
  0.00000000e+00,  0.00000000e+00,  1.13668486e+00,
  0.00000000e+00,  0.00000000e+00,  1.13668486e+00,
};

const double pylith::feassemble::IntegratorDataInertia2Din3DTwo::_valsLumped[] = {
  3.41005458e+00,  3.41005458e+00,  3.41005458e+00,
  6.82010916e+00,  6.82010916e+00,  6.82010916e+00,
  3.41005458e+00,  3.41005458e+00,  3.41005458e+00,
  3.41005458e+00,  3.41005458e+00,  3.41005458e+00,
  3.41005458e+00,  3.41005458e+00,  3.41005458e+00,
};

pylith::feassemble::IntegratorDataInertia2Din3DTwo::IntegratorDataInertia2Din3DTwo(void)
{ // constructor
  numVertices = _numVertices;
  spaceDim = _spaceDim;
  numCells = _numCells;
  cellDim = _cellDim;
  numCorners = _numCorners;
  numQuadPts = _numQuadPts;
  fiberDim = _fiberDim;
  vertices = const_cast<double*>(_vertices);
  cells = const_cast<int*>(_cells);
  quadPts = const_cast<double*>(_quadPts);
  quadWts = const_cast<double*>(_quadWts);
  basis = const_cast<double*>(_basis);
  basisDeriv = const_cast<double*>(_basisDeriv);
  fieldIn = const_cast<double*>(_fieldIn);
  valsAction = const_cast<double*>(_valsAction);
  valsMatrix = const_cast<double*>(_valsMatrix);
  valsLumped = const_cast<double*>(_valsLumped);
} // constructor

pylith::feassemble::IntegratorDataInertia2Din3DTwo::~IntegratorDataInertia2Din3DTwo(void)
{}


// End of file
