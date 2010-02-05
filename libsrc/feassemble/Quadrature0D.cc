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
  const int cellDim = 0;
  const int spaceDim = 1;
  const int numQuadPts = 1;
  const int numBasis = 1;

  assert(_quadRefCell.cellDim() == cellDim);
  assert(_quadRefCell.spaceDim() == spaceDim);
  assert(_quadRefCell.numQuadPts() == numQuadPts);
  assert(_quadRefCell.numBasis() == numBasis);
  assert(coordinatesCell.size() == numBasis*spaceDim);

  const double_array& basisDerivRef = _quadRefCell.basisDerivRef();

  zero();

  _quadPts = coordinatesCell;

  _jacobian[0] = 1.0;
  _jacobianDet[0] = 1.0;
  _jacobianInv[0] = 1.0;
  _basisDeriv[0] = basisDerivRef[0];

  PetscLogFlops(0);
} // computeGeometry


// End of file 
