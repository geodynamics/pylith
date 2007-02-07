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

#include "IntegratorDynExplicit.hh" // implementation of class methods

#include "Quadrature.hh"
#include "spatialdata/spatialdb/SpatialDB.hh"

#include <assert.h> // USES assert()

// ----------------------------------------------------------------------
// Constructor
pylith::feassemble::IntegratorDynExplicit::IntegratorDynExplicit(void) :
  _quadrature(0),
  _cellVector(0),
  _cellMatrix(0)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor
pylith::feassemble::IntegratorDynExplicit::~IntegratorDynExplicit(void)
{ // destructor
  delete _quadrature; _quadrature = 0;
  delete[] _cellVector; _cellVector = 0;
  delete[] _cellMatrix; _cellMatrix = 0;
} // destructor
  
// ----------------------------------------------------------------------
// Copy constructor
pylith::feassemble::IntegratorDynExplicit::IntegratorDynExplicit(const IntegratorDynExplicit& i) :
  _cellVector(0),
  _cellMatrix(0)
{ // copy constructor
  if (0 != i._quadrature)
    _quadrature = i._quadrature->clone();
} // copy constructor

// ----------------------------------------------------------------------
// Set quadrature for integrating finite-element quantities.
void
pylith::feassemble::IntegratorDynExplicit::quadrature(const Quadrature* q)
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
pylith::feassemble::IntegratorDynExplicit::_initCellVector(void)
{ // _initCellVector
  assert(0 != _quadrature);
  const int size = _quadrature->spaceDim() * _quadrature->numCorners();
  if (0 == _cellVector)
    _cellVector = (size > 0) ? new real_section_type::value_type[size] : 0;
  for (int i=0; i < size; ++i)
    _cellVector[i] = 0.0;
} // _initCellVector

// ----------------------------------------------------------------------
// Zero out vector containing result of integration actions for cell.
void
pylith::feassemble::IntegratorDynExplicit::_resetCellVector(void)
{ // _resetCellVector
  assert(0 != _quadrature);
  const int size = _quadrature->spaceDim() * _quadrature->numCorners();
  for (int i=0; i < size; ++i)
    _cellVector[i] = 0.0;
} // _resetCellVector

// ----------------------------------------------------------------------
// Initialize matrix containing result of integration for cell.
void
pylith::feassemble::IntegratorDynExplicit::_initCellMatrix(void)
{ // _initCellMatrix
  assert(0 != _quadrature);
  const int size =
    _quadrature->spaceDim() * _quadrature->numCorners() *
    _quadrature->spaceDim() * _quadrature->numCorners();
  if (0 == _cellMatrix)
    _cellMatrix = (size > 0) ? new real_section_type::value_type[size] : 0;
  for (int i=0; i < size; ++i)
    _cellMatrix[i] = 0.0;
} // _initCellMatrix

// ----------------------------------------------------------------------
// Zero out matrix containing result of integration for cell.
void
pylith::feassemble::IntegratorDynExplicit::_resetCellMatrix(void)
{ // _resetCellMatrix
  assert(0 != _quadrature);
  const int size =
    _quadrature->spaceDim() * _quadrature->numCorners() *
    _quadrature->spaceDim() * _quadrature->numCorners();
  for (int i=0; i < size; ++i)
    _cellMatrix[i] = 0.0;
} // _resetCellMatrix


// End of file 
