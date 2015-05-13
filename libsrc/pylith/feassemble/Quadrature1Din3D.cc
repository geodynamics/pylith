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

#include "Quadrature1Din3D.hh" // implementation of class methods

#include "QuadratureRefCell.hh" // USES QuadratureRefCell
#include "CellGeometry.hh" // USES CellGeometry

#include "petsc.h" // USES PetscLogFlops

#include <cassert> // USES assert()

#define ISOPARAMETRIC

// ----------------------------------------------------------------------
// Constructor
pylith::feassemble::Quadrature1Din3D::Quadrature1Din3D(const QuadratureRefCell& q) :
  QuadratureEngine(q)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor
pylith::feassemble::Quadrature1Din3D::~Quadrature1Din3D(void)
{ // destructor
} // destructor
  
// ----------------------------------------------------------------------
// Copy constructor.
pylith::feassemble::Quadrature1Din3D::Quadrature1Din3D(const Quadrature1Din3D& q) :
  QuadratureEngine(q)
{ // copy constructor
} // copy constructor

// ----------------------------------------------------------------------
// Compute geometric quantities for a cell at quadrature points.
void
pylith::feassemble::Quadrature1Din3D::computeGeometry(const PylithScalar* coordinatesCell,
						      const int coordinatesSize,
						      const int cell)
{ // computeGeometry
  assert(coordinatesCell);

  const int cellDim = 1;
  const int spaceDim = 3;

  const int numQuadPts = _quadRefCell.numQuadPts();
  const int numBasis = _quadRefCell.numBasis();

  const scalar_array& basis = _quadRefCell.basis();
  const scalar_array& basisDerivRef = _quadRefCell.basisDerivRef();

  assert(_quadRefCell.cellDim() == cellDim);
  assert(_quadRefCell.spaceDim() == spaceDim);
  assert(numBasis*spaceDim == coordinatesSize);

  zero();

  // Loop over quadrature points
  for (int iQuadPt=0; iQuadPt < numQuadPts; ++iQuadPt) {
    const int iQ = iQuadPt*numBasis;
    
    // Compute coordinates of quadrature point in cell
#if defined(ISOPARAMETRIC)
    // x = sum[i=0,n-1] (Ni * xi)
    // y = sum[i=0,n-1] (Ni * yi)
    // z = sum[i=0,n-1] (Ni * zi)
    for (int iBasis=0; iBasis < numBasis; ++iBasis) {
      const PylithScalar valueBasis = basis[iQ+iBasis];
      for (int iDim=0; iDim < spaceDim; ++iDim)
	_quadPts[iQuadPt*spaceDim+iDim] += 
	  valueBasis * coordinatesCell[iBasis*spaceDim+iDim];
    } // for
#else
    geometry.ptsRefToGlobal(&_quadPts[iQuadPt*spaceDim],
			    &quadPtsRef[iQuadPt*cellDim],
			    coordinatesCell, spaceDim, 1);
#endif
    
#if defined(ISOPARAMETRIC)
    // Compute Jacobian at quadrature point
    // J = [dx/dp
    //      dy/dp
    //      dz/dp]
    // dx/dp = sum[i=0,n-1] (dNi/dp * xi)
    // dy/dp = sum[i=0,n-1] (dNi/dp * yi)
    // dz/dp = sum[i=0,n-1] (dNi/dp * zi)
    for (int iBasis=0; iBasis < numBasis; ++iBasis) {
      const PylithScalar deriv = basisDerivRef[iQ+iBasis];
      for (int iDim=0; iDim < spaceDim; ++iDim)
	_jacobian[iQuadPt*spaceDim+iDim] += 
	  deriv * coordinatesCell[iBasis*spaceDim+iDim];
    } // for

    // Compute determinant of Jacobian at quadrature point
    // |J| = sqrt(transpose(J) J)
    PylithScalar det = 0.0;
    for (int iDim=0; iDim < spaceDim; ++iDim)
      det += _jacobian[iQuadPt*spaceDim+iDim] * 
	_jacobian[iQuadPt*spaceDim+iDim];
    det = sqrt(det);
    _checkJacobianDet(det, cell);
    _jacobianDet[iQuadPt] = det;
#else
    // Compute Jacobian and determinant of Jacobian at quadrature point
    geometry.jacobian(&_jacobian[iQuadPt*cellDim*spaceDim],
		      &_jacobianDet[iQuadPt],
		      coordinatesCell, &quadPtsRef[iQuadPt*cellDim],
		      spaceDim, 1);
    _checkJacobianDet(_jacobianDet[iQuadPt], cell);
#endif

    // Compute inverse of Jacobian at quadrature point
    // Jinv = 1.0/[J]
    for (int iDim=0; iDim < spaceDim; ++iDim)
      _jacobianInv[iQuadPt*spaceDim+iDim] = 
	1.0 / _jacobian[iQuadPt*spaceDim+iDim];

    // Compute derivatives of basis functions with respect to global
    // coordinates
    // dNi/dx = dNi/dp dp/dx + dNi/dq dq/dx + dNi/dr dr/dx
    const int iJ = iQuadPt*cellDim*spaceDim;
    for (int iBasis=0; iBasis < numBasis; ++iBasis) {
      const int iD = iQuadPt*numBasis*spaceDim + iBasis*spaceDim;
      const int iDR = iQuadPt*numBasis*cellDim + iBasis*cellDim;
      for (int iDim=0; iDim < spaceDim; ++iDim)
        for (int jDim=0; jDim < cellDim; ++jDim)
          _basisDeriv[iD+iDim] +=
	    basisDerivRef[iDR+jDim] * _jacobianInv[iJ+jDim*spaceDim+iDim];
    } // for
  } // for
  
  PetscLogFlops(numQuadPts * (1 + numBasis*spaceDim*2 +
			      spaceDim*1 +
			      numBasis*spaceDim*cellDim*2));

} // computeGeometry


// End of file 
