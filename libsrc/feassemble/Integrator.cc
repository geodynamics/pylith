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

#include "Quadrature.hh" // USES Quadrature
#include "pylith/utils/constdefs.h" // USES MAXDOUBLE

#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional
#include "spatialdata/spatialdb/GravityField.hh" // USES GravityField

#include <cassert> // USES assert()

// ----------------------------------------------------------------------
// Constructor
template<typename quadrature_type>
pylith::feassemble::Integrator<quadrature_type>::Integrator(void) :
  _dt(-1.0),
  _quadrature(0),
  _normalizer(new spatialdata::units::Nondimensional),
  _gravityField(0),
  _needNewJacobian(true),
  _useSolnIncr(false)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor
template<typename quadrature_type>
pylith::feassemble::Integrator<quadrature_type>::~Integrator(void)
{ // destructor
  delete _quadrature; _quadrature = 0;
  delete _normalizer; _normalizer = 0;
} // destructor
  
// ----------------------------------------------------------------------
// Set quadrature for integrating finite-element quantities.
template<typename quadrature_type>
void
pylith::feassemble::Integrator<quadrature_type>::quadrature(const quadrature_type* q)
{ // quadrature
  delete _quadrature;
  _quadrature = (0 != q) ? new quadrature_type(*q) : 0;

  // Deallocate cell vector and matrix since size may change
  _cellVector.resize(0);
  _cellMatrix.resize(0);
} // quadrature

// ----------------------------------------------------------------------
// Set manager of scales used to nondimensionalize problem.
template<typename quadrature_type>
void
pylith::feassemble::Integrator<quadrature_type>::normalizer(const spatialdata::units::Nondimensional& dim)
{ // normalizer
  if (0 == _normalizer)
    _normalizer = new spatialdata::units::Nondimensional(dim);
  else
    *_normalizer = dim;
} // normalizer

// ----------------------------------------------------------------------
// Set gravity field.
template<typename quadrature_type>
void
pylith::feassemble::Integrator<quadrature_type>::gravityField(spatialdata::spatialdb::GravityField* const gravityField)
{ // gravityField
  _gravityField = gravityField;
} // gravityField

// ----------------------------------------------------------------------
// Get stable time step for advancing from time t to time t+dt.
template<typename quadrature_type>
double
pylith::feassemble::Integrator<quadrature_type>::stableTimeStep(void) const
{ // stableTimeStep
  // Assume any time step will work.
  return pylith::PYLITH_MAXDOUBLE;
} // stableTimeStep

// ----------------------------------------------------------------------
// Initialize vector containing result of integration action for cell.
template<typename quadrature_type>
void
pylith::feassemble::Integrator<quadrature_type>::_initCellVector(void)
{ // _initCellVector
  assert(0 != _quadrature);
  const int size = _quadrature->spaceDim() * _quadrature->numBasis();
  _cellVector.resize(size);
  _cellVector = 0.0;
} // _initCellVector

// ----------------------------------------------------------------------
// Zero out vector containing result of integration actions for cell.
template<typename quadrature_type>
void
pylith::feassemble::Integrator<quadrature_type>::_resetCellVector(void)
{ // _resetCellVector
  _cellVector = 0.0;
} // _resetCellVector

// ----------------------------------------------------------------------
// Initialize matrix containing result of integration for cell.
template<typename quadrature_type>
void
pylith::feassemble::Integrator<quadrature_type>::_initCellMatrix(void)
{ // _initCellMatrix
  assert(0 != _quadrature);
  const int size =
    _quadrature->spaceDim() * _quadrature->numBasis() *
    _quadrature->spaceDim() * _quadrature->numBasis();
  _cellMatrix.resize(size);
  _cellMatrix = 0.0;
} // _initCellMatrix

// ----------------------------------------------------------------------
// Zero out matrix containing result of integration for cell.
template<typename quadrature_type>
void
pylith::feassemble::Integrator<quadrature_type>::_resetCellMatrix(void)
{ // _resetCellMatrix
  _cellMatrix = 0.0;
} // _resetCellMatrix


// End of file 
