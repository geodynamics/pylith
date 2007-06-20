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

#include <portinfo>

#include "Integrator.hh" // implementation of class methods

#include "Quadrature.hh" // USES Quadrature

#include <assert.h> // USES assert()

// ----------------------------------------------------------------------
// Constructor
pylith::feassemble::Integrator::Integrator(void) :
  _dt(-1.0),
  _quadrature(0),
  _cellVector(0),
  _cellMatrix(0),
  _needNewJacobian(true),
  _useSolnIncr(false)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor
pylith::feassemble::Integrator::~Integrator(void)
{ // destructor
  delete _quadrature; _quadrature = 0;
  delete[] _cellVector; _cellVector = 0;
  delete[] _cellMatrix; _cellMatrix = 0;
} // destructor
  
// ----------------------------------------------------------------------
// Set quadrature for integrating finite-element quantities.
void
pylith::feassemble::Integrator::quadrature(const Quadrature* q)
{ // quadrature
  delete _quadrature;
  _quadrature = (0 != q) ? q->clone() : 0;

  // Deallocate cell vector and matrix since size may change
  delete[] _cellVector; _cellVector = 0;
  delete[] _cellMatrix; _cellMatrix = 0;
} // quadrature

// ----------------------------------------------------------------------
// Initialize vector containing result of integration action for cell.
void
pylith::feassemble::Integrator::_initCellVector(void)
{ // _initCellVector
  assert(0 != _quadrature);
  const int size = _quadrature->spaceDim() * _quadrature->numBasis();
  if (0 == _cellVector)
    _cellVector = (size > 0) ? new real_section_type::value_type[size] : 0;
  for (int i=0; i < size; ++i)
    _cellVector[i] = 0.0;
} // _initCellVector

// ----------------------------------------------------------------------
// Zero out vector containing result of integration actions for cell.
void
pylith::feassemble::Integrator::_resetCellVector(void)
{ // _resetCellVector
  assert(0 != _quadrature);
  assert(0 != _cellVector);
  const int size = _quadrature->spaceDim() * _quadrature->numBasis();
  for (int i=0; i < size; ++i)
    _cellVector[i] = 0.0;
} // _resetCellVector

// ----------------------------------------------------------------------
// Initialize matrix containing result of integration for cell.
void
pylith::feassemble::Integrator::_initCellMatrix(void)
{ // _initCellMatrix
  assert(0 != _quadrature);
  const int size =
    _quadrature->spaceDim() * _quadrature->numBasis() *
    _quadrature->spaceDim() * _quadrature->numBasis();
  if (0 == _cellMatrix)
    _cellMatrix = (size > 0) ? new real_section_type::value_type[size] : 0;
  for (int i=0; i < size; ++i)
    _cellMatrix[i] = 0.0;
} // _initCellMatrix

// ----------------------------------------------------------------------
// Zero out matrix containing result of integration for cell.
void
pylith::feassemble::Integrator::_resetCellMatrix(void)
{ // _resetCellMatrix
  assert(0 != _quadrature);
  assert(0 != _cellMatrix);
  const int size =
    _quadrature->spaceDim() * _quadrature->numBasis() *
    _quadrature->spaceDim() * _quadrature->numBasis();
  for (int i=0; i < size; ++i)
    _cellMatrix[i] = 0.0;
} // _resetCellMatrix


// End of file 
