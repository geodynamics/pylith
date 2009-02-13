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

#include "Quadrature1Din2D.hh" // implementation of class methods

#include "QuadratureRefCell.hh" // USES QuadratureRefCell
#include "CellGeometry.hh" // USES CellGeometry

#include "petsc.h" // USES PetscLogFlops

#include <cassert> // USES assert()

#define ISOPARAMETRIC

// ----------------------------------------------------------------------
// Constructor
pylith::feassemble::Quadrature1Din2D::Quadrature1Din2D(const QuadratureRefCell& q) :
  QuadratureEngine(q)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor
pylith::feassemble::Quadrature1Din2D::~Quadrature1Din2D(void)
{ // destructor
} // destructor
  
// ----------------------------------------------------------------------
// Copy constructor.
pylith::feassemble::Quadrature1Din2D::Quadrature1Din2D(const Quadrature1Din2D& q) :
  QuadratureEngine(q)
{ // copy constructor
} // copy constructor

// ----------------------------------------------------------------------
// Compute geometric quantities for a cell at quadrature points.
void
pylith::feassemble::Quadrature1Din2D::computeGeometry(const double* vertCoords,
						      const int coordDim,
						      const int cell)
{ // computeGeometry
  assert(0 != vertCoords);
  assert(2 == coordDim);

  const int cellDim = _quadRefCell.cellDim();
  const int spaceDim = _quadRefCell.spaceDim();
  const int numQuadPts = _quadRefCell.numQuadPts();
  const int numBasis = _quadRefCell.numBasis();

  const double_array& basis = _quadRefCell.basis();
  const double_array& quadPtsRef = _quadRefCell.quadPtsRef();
  const double_array& basisDerivRef = _quadRefCell.basisDerivRef();
  const CellGeometry& geometry = _quadRefCell.refGeometry();

  assert(1 == cellDim);
  assert(2 == spaceDim);
  zero();

  // Loop over quadrature points
  for (int iQuadPt=0; iQuadPt < numQuadPts; ++iQuadPt) {
    
    // Compute coordinates of quadrature point in cell
#if defined(ISOPARAMETRIC)
    // x = sum[i=0,n-1] (Ni * xi)
    // y = sum[i=0,n-1] (Ni * yi)
    for (int iBasis=0; iBasis < numBasis; ++iBasis) {
      const double valueBasis = basis[iQuadPt*numBasis+iBasis];
      for (int iDim=0; iDim < spaceDim; ++iDim)
	_quadPts[iQuadPt*spaceDim+iDim] +=
	  valueBasis * vertCoords[iBasis*spaceDim+iDim];
    } // for
#else
    geometry.coordsRefToGlobal(&_quadPts[iQuadPt*spaceDim],
			       &quadPtsRef[iQuadPt*cellDim],
			       vertCoords, spaceDim);
#endif

#if defined(ISOPARAMETRIC)
    // Compute Jacobian at quadrature point
    // J = [dx/dp
    //      dy/dp]
    // dx/dp = sum[i=0,n-1] (dNi/dp * xi)
    // dy/dp = sum[i=0,n-1] (dNi/dp * yi)
    for (int iBasis=0; iBasis < numBasis; ++iBasis) {
      const double deriv = basisDerivRef[iQuadPt*numBasis+iBasis];
      for (int iDim=0; iDim < spaceDim; ++iDim)
	_jacobian[iQuadPt*spaceDim+iDim] += 
	  deriv * vertCoords[iBasis*spaceDim+iDim];
    } // for

    // Compute determinant of Jacobian at quadrature point
    // |J| = sqrt(transpose(J) J)
    double det = 0.0;
    for (int iDim=0; iDim < spaceDim; ++iDim)
      det += _jacobian[iQuadPt*spaceDim+iDim] * 
	_jacobian[iQuadPt*spaceDim+iDim];
    det = sqrt(det);
    _checkJacobianDet(det, cell);
    _jacobianDet[iQuadPt] = det;
#else
    // Compute Jacobian and determinant of Jacobian at quadrature point
    geometry.jacobian(&_jacobian[iQuadPt*_cellDim*spaceDim],
		      &_jacobianDet[iQuadPt],
		      vertCoords, &quadPtsRef[iQuadPt*cellDim], spaceDim);
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
    for (int iBasis=0; iBasis < numBasis; ++iBasis)
      for (int iDim=0; iDim < spaceDim; ++iDim)
	for (int jDim=0; jDim < cellDim; ++jDim)
	  _basisDeriv[iQuadPt*numBasis*spaceDim+iBasis*spaceDim+iDim] +=
	    basisDerivRef[iQuadPt*numBasis*cellDim+iBasis*cellDim+jDim] *
	    _jacobianInv[iQuadPt*cellDim*spaceDim+jDim*spaceDim+iDim];
  } // for

  PetscLogFlops(numQuadPts * (1 + numBasis*spaceDim*2 +
			      spaceDim*1 +
			      numBasis*spaceDim*cellDim*2));
} // computeGeometry


// End of file 
