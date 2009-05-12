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

#include "Quadrature0D.hh" // implementation of class methods

#include "QuadratureRefCell.hh" // USES QuadratureRefCell
#include "CellGeometry.hh" // USES CellGeometry

#include "petsc.h" // USES PetscLogFlops()

#include <cassert> // USES assert()

// ----------------------------------------------------------------------
// Constructor
pylith::feassemble::Quadrature0D::Quadrature0D(const QuadratureRefCell& q) :
  QuadratureEngine(q)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor
pylith::feassemble::Quadrature0D::~Quadrature0D(void)
{ // destructor
} // destructor
  
// ----------------------------------------------------------------------
// Copy constructor.
pylith::feassemble::Quadrature0D::Quadrature0D(const Quadrature0D& q) :
  QuadratureEngine(q)
{ // copy constructor
} // copy constructor

// ----------------------------------------------------------------------
// Compute geometric quantities for a cell at quadrature points.
void
pylith::feassemble::Quadrature0D::computeGeometry(const double_array& coordinatesCell,
						  const int cell)
{ // computeGeometry
  const int cellDim = _quadRefCell.cellDim();
  const int spaceDim = _quadRefCell.spaceDim();
  const int numQuadPts = _quadRefCell.numQuadPts();
  const int numBasis = _quadRefCell.numBasis();
  assert(coordinatesCell.size() == numBasis*spaceDim);

  const double_array& basisDerivRef = _quadRefCell.basisDerivRef();

  assert(0 == cellDim);
  assert(1 == numQuadPts);
  assert(1 == numBasis);

  zero();
  assert(numBasis*1 == coordinatesCell.size());

  _quadPts = coordinatesCell;

  _jacobian[0] = 1.0;
  _jacobianDet[0] = 1.0;
  _jacobianInv[0] = 1.0;
  _basisDeriv[0] = basisDerivRef[0];

  PetscLogFlops(0);
} // computeGeometry


// End of file 
