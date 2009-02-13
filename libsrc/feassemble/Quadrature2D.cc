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

#include "Quadrature2D.hh" // implementation of class methods

#include "QuadratureRefCell.hh" // USES QuadratureRefCell
#include "CellGeometry.hh" // USES CellGeometry

#include "petsc.h" // USES PetscLogFlops

#include <cassert> // USES assert()

#define ISOPARAMETRIC

// ----------------------------------------------------------------------
// Constructor
pylith::feassemble::Quadrature2D::Quadrature2D(const QuadratureRefCell& q) :
  QuadratureEngine(q)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor
pylith::feassemble::Quadrature2D::~Quadrature2D(void)
{ // destructor
} // destructor
  
// ----------------------------------------------------------------------
// Copy constructor.
pylith::feassemble::Quadrature2D::Quadrature2D(const Quadrature2D& q) :
  QuadratureEngine(q)
{ // copy constructor
} // copy constructor

// ----------------------------------------------------------------------
// Compute geometric quantities for a cell at quadrature points.
void
pylith::feassemble::Quadrature2D::computeGeometry(const double* vertCoords,
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

  assert(2 == cellDim);
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
    // J = [dx/dp dx/dq]
    //     [dy/dp dy/dq]
    // dx/dp = sum[i=0,n-1] (dNi/dp * xi)
    // dx/dq = sum[i=0,n-1] (dNi/dq * xi)
    // dy/dp = sum[i=0,n-1] (dNi/dp * yi)
    // dy/dq = sum[i=0,n-1] (dNi/dq * yi)
    for (int iBasis=0; iBasis < numBasis; ++iBasis)
      for (int iCol=0; iCol < cellDim; ++iCol) {
	const double deriv = 
	  basisDerivRef[iQuadPt*numBasis*spaceDim+iBasis*cellDim+iCol];
	for (int iRow=0; iRow < spaceDim; ++iRow)
	  _jacobian[iQuadPt*cellDim*spaceDim+iRow*cellDim+iCol] +=
	    deriv * vertCoords[iBasis*spaceDim+iRow];
      } // for
  
    // Compute determinant of Jacobian at quadrature point
    // |J| = j00*j11-j01*j10
    const int iJ = iQuadPt*cellDim*spaceDim;
    const int i00 = iJ + 0*spaceDim + 0;
    const int i01 = iJ + 0*spaceDim + 1;
    const int i10 = iJ + 1*spaceDim + 0;
    const int i11 = iJ + 1*spaceDim + 1;
    const double det = 
      _jacobian[i00]*_jacobian[i11] - 
      _jacobian[i01]*_jacobian[i10];
    _checkJacobianDet(det, cell);
    _jacobianDet[iQuadPt] = det;
#else
    // Compute Jacobian and determinant of Jacobian at quadrature point
    geometry.jacobian(&_jacobian[iQuadPt*cellDim*spaceDim],
		      &_jacobianDet[iQuadPt],
		      vertCoords, &quadPtsRef[iQuadPt*cellDim], spaceDim);
    _checkJacobianDet(_jacobianDet[iQuadPt], cell);

    const int iJ = iQuadPt*cellDim*spaceDim;
    const int i00 = iJ + 0*spaceDim + 0;
    const int i01 = iJ + 0*spaceDim + 1;
    const int i10 = iJ + 1*spaceDim + 0;
    const int i11 = iJ + 1*spaceDim + 1;
    const double det = _jacobianDet[iQuadPt];
#endif

    // Compute inverse of Jacobian at quadrature point
    // Jinv = 1/det*[ j11 -j01]
    //              [-j10  j00]
    _jacobianInv[i00] =  _jacobian[i11] / det;
    _jacobianInv[i01] = -_jacobian[i01] / det;
    _jacobianInv[i10] = -_jacobian[i10] / det;
    _jacobianInv[i11] =  _jacobian[i00] / det;

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

  PetscLogFlops(numQuadPts*(4 +
			    numBasis*spaceDim*2 +
			    numBasis*spaceDim*cellDim*2));
} // computeGeometry


// End of file 
