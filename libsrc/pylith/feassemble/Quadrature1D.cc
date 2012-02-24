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
// Copyright (c) 2010-2012 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

#include <portinfo>

#include "Quadrature1D.hh" // implementation of class methods

#include "QuadratureRefCell.hh" // USES QuadratureRefCell
#include "CellGeometry.hh" // USES CellGeometry

#include "petsc.h" // USES PetscLogFlops

#include <cassert> // USES assert()

#define ISOPARAMETRIC

// ----------------------------------------------------------------------
// Constructor
pylith::feassemble::Quadrature1D::Quadrature1D(const QuadratureRefCell& q) :
  QuadratureEngine(q)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor
pylith::feassemble::Quadrature1D::~Quadrature1D(void)
{ // destructor
} // destructor
  
// ----------------------------------------------------------------------
// Copy constructor.
pylith::feassemble::Quadrature1D::Quadrature1D(const Quadrature1D& q) :
  QuadratureEngine(q)
{ // copy constructor
} // copy constructor

// ----------------------------------------------------------------------
// Compute geometric quantities for a cell at quadrature points.
void
pylith::feassemble::Quadrature1D::computeGeometry(const scalar_array& coordinatesCell,
						  const int cell)
{ // computeGeometry
  const int cellDim = 1;
  const int spaceDim = 1;

  const int numQuadPts = _quadRefCell.numQuadPts();
  const int numBasis = _quadRefCell.numBasis();

  const scalar_array& basis = _quadRefCell.basis();
  const scalar_array& quadPtsRef = _quadRefCell.quadPtsRef();
  const scalar_array& basisDerivRef = _quadRefCell.basisDerivRef();
  const CellGeometry& geometry = _quadRefCell.refGeometry();

  assert(_quadRefCell.cellDim() == cellDim);
  assert(_quadRefCell.spaceDim() == spaceDim);
  assert(numBasis*spaceDim == coordinatesCell.size());

  zero();
  
  // Loop over quadrature points
  for (int iQuadPt=0; iQuadPt < numQuadPts; ++iQuadPt) {
    const int iQ = iQuadPt*numBasis;

    // Compute coordinates of quadrature point in cell
#if defined(ISOPARAMETRIC)
    // x = sum[i=0,n-1] (Ni * xi)
    for (int iBasis=0; iBasis < numBasis; ++iBasis)
      _quadPts[iQuadPt] += basis[iQ+iBasis]*coordinatesCell[iBasis];
#else
    geometry.ptsRefToGlobal(&_quadPts[iQuadPt], &quadPtsRef[iQuadPt],
			    &coordinatesCell[0], spaceDim, 1);
#endif

#if defined(ISOPARAMETRIC)
    // Compute Jacobian at quadrature point
    // J = dx/dp = sum[i=0,n-1] (dNi/dp * xi)
    for (int iBasis=0; iBasis < numBasis; ++iBasis)
      _jacobian[iQuadPt] += basisDerivRef[iQ+iBasis] * coordinatesCell[iBasis];

    // Compute determinant of Jacobian at quadrature point
    // |J| = j00
    const PylithScalar det = _jacobian[iQuadPt];
    _checkJacobianDet(det, cell);
    _jacobianDet[iQuadPt] = _jacobian[iQuadPt];
#else
    // Compute Jacobian and determinant of Jacobian at quadrature point
    geometry.jacobian(&_jacobian[iQuadPt], &_jacobianDet[iQuadPt],
		      &coordinatesCell[0], &quadPtsRef[iQuadPt], spaceDim, 1);
    _checkJacobianDet(_jacobianDet[iQuadPt], cell);
#endif
    
    // Compute inverse of Jacobian at quadrature point
    // Jinv = 1/j00
    _jacobianInv[iQuadPt] = 1.0 / _jacobianDet[iQuadPt];

    assert(numQuadPts*numBasis*spaceDim == _basisDeriv.size());
    assert(numQuadPts*numBasis*cellDim == basisDerivRef.size());
    assert(numQuadPts*cellDim*spaceDim == _jacobianInv.size());

    // Compute derivatives of basis functions with respect to global
    // coordinates
    // dNi/dx = dNi/dp dp/dx + dNi/dq dq/dx + dNi/dr dr/dx
    for (int iBasis=0; iBasis < numBasis; ++iBasis)
      _basisDeriv[iQ+iBasis] +=
          basisDerivRef[iQ+iBasis] * _jacobianInv[iQuadPt];
  } // for

  PetscLogFlops(numQuadPts * (1 + numBasis * 4));
} // computeGeometry


// End of file 
