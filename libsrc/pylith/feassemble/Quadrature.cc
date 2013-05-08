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
// Copyright (c) 2010-2013 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

#include <portinfo>

#include "CellGeometry.hh" // USES CellGeometry

#include "QuadratureEngine.hh" // USES QuadratureEngine
#include "Quadrature0D.hh"
#include "Quadrature1D.hh"
#include "Quadrature1Din2D.hh"
#include "Quadrature1Din3D.hh"
#include "Quadrature2D.hh"
#include "Quadrature2Din3D.hh"
#include "Quadrature3D.hh"

#include <cassert> // USES assert()
#include <stdexcept> // USES std::runtime_error
#include <sstream> // USES std::ostringstream

// ----------------------------------------------------------------------
// Constructor
template<typename mesh_type>
pylith::feassemble::Quadrature<mesh_type>::Quadrature(void) :
  _engine(0),
  _checkConditioning(false)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor
template<typename mesh_type>
pylith::feassemble::Quadrature<mesh_type>::~Quadrature(void)
{ // destructor
  deallocate();
} // destructor
  
// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
template<typename mesh_type>
void
pylith::feassemble::Quadrature<mesh_type>::deallocate(void)
{ // deallocate
  PYLITH_METHOD_BEGIN;

  QuadratureRefCell::deallocate();

  delete _engine; _engine = 0;

  PYLITH_METHOD_END;
} // deallocate
  
// ----------------------------------------------------------------------
// Copy constructor
template<typename mesh_type>
pylith::feassemble::Quadrature<mesh_type>::Quadrature(const Quadrature& q) :
  QuadratureRefCell(q),
  _engine(0),
  _checkConditioning(q._checkConditioning)
{ // copy constructor
  PYLITH_METHOD_BEGIN;

  if (q._engine)
    _engine = q._engine->clone();

  PYLITH_METHOD_END;
} // copy constructor

// ----------------------------------------------------------------------
// Setup quadrature engine.
template<typename mesh_type>
void
pylith::feassemble::Quadrature<mesh_type>::initializeGeometry(void)
{ // initializeGeometry
  PYLITH_METHOD_BEGIN;

  clear();
  assert(!_engine);

  const int cellDim = _cellDim;
  const int spaceDim = _spaceDim;

  if (1 == spaceDim)
    if (1 == cellDim)
      _engine = new Quadrature1D(*this);
    else if (0 == cellDim)
      _engine = new Quadrature0D(*this);
    else {
      std::cerr << "Unknown quadrature case with cellDim '" 
		<< cellDim << "' and spaceDim '" << spaceDim << "'" 
		<< std::endl;
      assert(0);
    } // if/else
  else if (2 == spaceDim)
    if (2 == cellDim)
      _engine = new Quadrature2D(*this);
    else if (1 == cellDim)
      _engine = new Quadrature1Din2D(*this);
    else if (0 == cellDim)
      _engine = new Quadrature0D(*this);
    else {
      std::cerr << "Unknown quadrature case with cellDim '" 
		<< cellDim << "' and spaceDim '" << spaceDim << "'" 
		<< std::endl;
      assert(0);
    } // if/else
  else if (3 == spaceDim)
    if (3 == cellDim)
      _engine = new Quadrature3D(*this);
    else if (2 == cellDim)
      _engine = new Quadrature2Din3D(*this);
    else if (1 == cellDim)
      _engine = new Quadrature1Din3D(*this);
    else if (0 == cellDim)
      _engine = new Quadrature0D(*this);
    else {
      std::cerr << "Unknown quadrature case with cellDim '" 
		<< cellDim << "' and spaceDim '" << spaceDim << "'" 
		<< std::endl;
      assert(0);
    } // if/else
  else {
    std::cerr << "Unknown quadrature case with cellDim '" 
	      << cellDim << "' and spaceDim '" << spaceDim << "'" 
	      << std::endl;
    assert(0);
  } // if/else

  assert(_engine);
  _engine->initialize();

  PYLITH_METHOD_END;
} // initializeGeometry

// ----------------------------------------------------------------------
// Deallocate temporary storage;
template<typename mesh_type>
void
pylith::feassemble::Quadrature<mesh_type>::clear(void)
{ // clear
  PYLITH_METHOD_BEGIN;

  delete _engine; _engine = 0;

  PYLITH_METHOD_END;
} // clear


// End of file 
