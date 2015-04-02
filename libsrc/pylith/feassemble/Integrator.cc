// -*- C++ -*-
//
// ======================================================================
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2015 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

#include <portinfo>

#include "Integrator.hh" // Implementation of class methods

#include "Quadrature.hh" // USES Quadrature

#include "pylith/utils/error.h" // USES PYLITH_METHOD_BEGIN/END
#include "pylith/utils/EventLogger.hh" // USES EventLogger
#include "pylith/utils/constdefs.h" // USES MAXSCALAR

#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional
#include "spatialdata/spatialdb/GravityField.hh" // USES GravityField

#include <cassert> // USES assert()

// ----------------------------------------------------------------------
// Constructor
pylith::feassemble::Integrator::Integrator(void) :
  _dt(-1.0),
  _quadrature(0),
  _normalizer(new spatialdata::units::Nondimensional),
  _gravityField(0),
  _logger(0),
  _needNewJacobian(true),
  _isJacobianSymmetric(true)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor
pylith::feassemble::Integrator::~Integrator(void)
{ // destructor
  deallocate();
} // destructor
  
// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::feassemble::Integrator::deallocate(void)
{ // deallocate
  PYLITH_METHOD_BEGIN;

  delete _quadrature; _quadrature = 0;
  delete _normalizer; _normalizer = 0;
  delete _logger; _logger = 0;
  _gravityField = 0; /// Memory managed elsewhere :TODO: use shared pointer

  PYLITH_METHOD_END;
} // deallocate
  
// ----------------------------------------------------------------------
// Get quadrature for integrating finite-element quantities.
const pylith::feassemble::Quadrature&
pylith::feassemble::Integrator::quadrature()
{ // quadrature
  return *_quadrature;
} // quadrature
  
// ----------------------------------------------------------------------
// Set quadrature for integrating finite-element quantities.
void
pylith::feassemble::Integrator::quadrature(const Quadrature* q)
{ // quadrature
  delete _quadrature; _quadrature = (q) ? new Quadrature(*q) : 0;

  // Deallocate cell vector and matrix since size may change
  _cellVector.resize(0);
  _cellMatrix.resize(0);
} // quadrature

// ----------------------------------------------------------------------
// Set manager of scales used to nondimensionalize problem.
void
pylith::feassemble::Integrator::normalizer(const spatialdata::units::Nondimensional& dim)
{ // normalizer
  if (0 == _normalizer)
    _normalizer = new spatialdata::units::Nondimensional(dim);
  else
    *_normalizer = dim;
} // normalizer

// ----------------------------------------------------------------------
// Set gravity field.
void
pylith::feassemble::Integrator::gravityField(spatialdata::spatialdb::GravityField* const gravityField)
{ // gravityField
  _gravityField = gravityField;
} // gravityField

// ----------------------------------------------------------------------
// Get stable time step for advancing from time t to time t+dt.
PylithScalar
pylith::feassemble::Integrator::stableTimeStep(const topology::Mesh& mesh)
{ // stableTimeStep
  // Assume any time step will work.
  return pylith::PYLITH_MAXSCALAR;
} // stableTimeStep

// ----------------------------------------------------------------------
// Initialize vector containing result of integration action for cell.
void
pylith::feassemble::Integrator::_initCellVector(void)
{ // _initCellVector
  assert(_quadrature);
  const int size = _quadrature->spaceDim() * _quadrature->numBasis();
  _cellVector.resize(size);
  _cellVector = 0.0;
} // _initCellVector

// ----------------------------------------------------------------------
// Zero out vector containing result of integration actions for cell.
void
pylith::feassemble::Integrator::_resetCellVector(void)
{ // _resetCellVector
  _cellVector = 0.0;
} // _resetCellVector

// ----------------------------------------------------------------------
// Initialize matrix containing result of integration for cell.
void
pylith::feassemble::Integrator::_initCellMatrix(void)
{ // _initCellMatrix
  assert(_quadrature);
  const int size =
    _quadrature->spaceDim() * _quadrature->numBasis() *
    _quadrature->spaceDim() * _quadrature->numBasis();
  _cellMatrix.resize(size);
  _cellMatrix = 0.0;
} // _initCellMatrix

// ----------------------------------------------------------------------
// Zero out matrix containing result of integration for cell.
void
pylith::feassemble::Integrator::_resetCellMatrix(void)
{ // _resetCellMatrix
  _cellMatrix = 0.0;
} // _resetCellMatrix

// ----------------------------------------------------------------------
// Lump cell matrix, putting the result in the cell vector using
// equivalent forces for rigid body motion.
void
pylith::feassemble::Integrator::_lumpCellMatrix(void)
{ // _lumpCellMatrix
  assert(_quadrature);
  const int numBasis = _quadrature->numBasis();
  const int spaceDim = _quadrature->spaceDim();

  _cellVector = 0.0;
  PylithScalar value = 0.0;
  for (int iBasis=0; iBasis < numBasis; ++iBasis)
    for (int iDim=0; iDim < spaceDim; ++iDim) {
      value = 0.0;
      const int indexR = (iBasis*spaceDim+iDim) * numBasis * spaceDim;
      for (int jBasis=0; jBasis < numBasis; ++jBasis)
	value += _cellMatrix[indexR+jBasis*spaceDim+iDim];
      if (value < 0.0) {
	throw std::runtime_error("Negative diagonal entry computed when "
				 "lumping Jacobian matrix.");
      } // for
      _cellVector[iBasis*spaceDim+iDim] = value;
    } // for

  PetscLogFlops(numBasis*numBasis*spaceDim);
} // _lumpCellMatrix


// End of file 
