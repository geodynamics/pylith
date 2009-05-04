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

#include "Quadrature3D.hh" // implementation of class methods

#include "QuadratureRefCell.hh" // USES QuadratureRefCell
#include "CellGeometry.hh" // USES CellGeometry

#include "petsc.h" // USES PetscLogFlops

#include <cassert> // USES assert()

#define ISOPARAMETRIC

// ----------------------------------------------------------------------
// Constructor
pylith::feassemble::Quadrature3D::Quadrature3D(const QuadratureRefCell& q) :
  QuadratureEngine(q)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor
pylith::feassemble::Quadrature3D::~Quadrature3D(void)
{ // destructor
} // destructor
  
// ----------------------------------------------------------------------
// Copy constructor.
pylith::feassemble::Quadrature3D::Quadrature3D(const Quadrature3D& q) :
  QuadratureEngine(q)
{ // copy constructor
} // copy constructor

// ----------------------------------------------------------------------
// Compute geometric quantities for a cell at quadrature points.
void
pylith::feassemble::Quadrature3D::computeGeometry(const double* vertCoords,
						  const int coordDim,
						  const int cell)
{ // computeGeometry
  assert(0 != vertCoords);
  assert(3 == coordDim);

  const int cellDim = _quadRefCell.cellDim();
  const int spaceDim = _quadRefCell.spaceDim();
  const int numQuadPts = _quadRefCell.numQuadPts();
  const int numBasis = _quadRefCell.numBasis();

  const double_array& basis = _quadRefCell.basis();
  const double_array& quadPtsRef = _quadRefCell.quadPtsRef();
  const double_array& basisDerivRef = _quadRefCell.basisDerivRef();
  const CellGeometry& geometry = _quadRefCell.refGeometry();

  assert(3 == cellDim);
  assert(3 == spaceDim);
  zero();
  
  // Loop over quadrature points
  for (int iQuadPt=0; iQuadPt < numQuadPts; ++iQuadPt) {
    
    // Compute coordinates of quadrature point in cell
#if defined(ISOPARAMETRIC)
    // x = sum[i=0,n-1] (Ni * xi)
    // y = sum[i=0,n-1] (Ni * yi)
    // z = sum[i=0,n-1] (Ni * zi)
    for (int iBasis=0; iBasis < numBasis; ++iBasis) {
      const double valueBasis = basis[iQuadPt*numBasis+iBasis];
      for (int iDim=0; iDim < spaceDim; ++iDim)
	_quadPts[iQuadPt*spaceDim+iDim] += 
	  valueBasis * vertCoords[iBasis*spaceDim+iDim];
    } // for
#else
    geometry.coordsRefToGlobal(&quadPts[iQuadPt*spaceDim],
			       &quadPtsRef[iQuadPt*cellDim],
			       vertCoords, spaceDim);
#endif
    
#if defined(ISOPARAMETRIC)
    // Compute Jacobian at quadrature point
    // J = [dx/dp dx/dq dx/dr]
    //     [dy/dp dy/dq dy/dr]
    //     [dz/dp dz/dq dz/dr]
    // dx/dp = sum[i=0,n-1] (dNi/dp * xi)
    // dx/dq = sum[i=0,n-1] (dNi/dq * xi)
    // dx/dr = sum[i=0,n-1] (dNi/dr * xi)
    // dy/dp = sum[i=0,n-1] (dNi/dp * yi)
    // dy/dq = sum[i=0,n-1] (dNi/dq * yi)
    // dy/dr = sum[i=0,n-1] (dNi/dr * yi)
    // dz/dp = sum[i=0,n-1] (dNi/dp * zi)
    // dz/dq = sum[i=0,n-1] (dNi/dq * zi)
    // dz/dr = sum[i=0,n-1] (dNi/dr * zi)
    for (int iBasis=0; iBasis < numBasis; ++iBasis)
      for (int iCol=0; iCol < cellDim; ++iCol) {
	const double deriv = 
	  basisDerivRef[iQuadPt*numBasis*spaceDim+iBasis*cellDim+iCol];
	for (int iRow=0; iRow < spaceDim; ++iRow)
	  _jacobian[iQuadPt*cellDim*spaceDim+iRow*cellDim+iCol] += 
	    deriv * vertCoords[iBasis*spaceDim+iRow];
      } // for

    // Compute determinant of Jacobian at quadrature point
    // |J| = j00*(j11*j22-j12*j21) +
    //      -j01*(j10*j22-j12*j20) +
    //       j02*(j10*j21-j11*j20)
    const int iJ = iQuadPt*cellDim*spaceDim;
    const int i00 = iJ + 0*spaceDim + 0;
    const int i01 = iJ + 0*spaceDim + 1;
    const int i02 = iJ + 0*spaceDim + 2;
    const int i10 = iJ + 1*spaceDim + 0;
    const int i11 = iJ + 1*spaceDim + 1;
    const int i12 = iJ + 1*spaceDim + 2;
    const int i20 = iJ + 2*spaceDim + 0;
    const int i21 = iJ + 2*spaceDim + 1;
    const int i22 = iJ + 2*spaceDim + 2;
    const double det = 
      _jacobian[i00]*(_jacobian[i11]*_jacobian[i22] -
		      _jacobian[i12]*_jacobian[i21]) -
      _jacobian[i01]*(_jacobian[i10]*_jacobian[i22] -
		      _jacobian[i12]*_jacobian[i20]) +
      _jacobian[i02]*(_jacobian[i10]*_jacobian[i21] -
		      _jacobian[i11]*_jacobian[i20]);
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
    const int i02 = iJ + 0*spaceDim + 2;
    const int i10 = iJ + 1*spaceDim + 0;
    const int i11 = iJ + 1*spaceDim + 1;
    const int i12 = iJ + 1*spaceDim + 2;
    const int i20 = iJ + 2*spaceDim + 0;
    const int i21 = iJ + 2*spaceDim + 1;
    const int i22 = iJ + 2*spaceDim + 2;
    const double det = _jacobianDet[iQuadPt];
#endif
    
    // Compute inverse of Jacobian at quadrature point
    _jacobianInv[i00] = (_jacobian[i11]*_jacobian[i22] -
			 _jacobian[i12]*_jacobian[i21]) / det;
    _jacobianInv[i01] = (_jacobian[i02]*_jacobian[i21] -
			 _jacobian[i01]*_jacobian[i22]) / det;
    _jacobianInv[i02] = (_jacobian[i01]*_jacobian[i12] -
			 _jacobian[i02]*_jacobian[i11]) / det;
    _jacobianInv[i10] = (_jacobian[i12]*_jacobian[i20] -
			 _jacobian[i10]*_jacobian[i22]) / det;
    _jacobianInv[i11] = (_jacobian[i00]*_jacobian[i22] -
			 _jacobian[i02]*_jacobian[i20]) / det;
    _jacobianInv[i12] = (_jacobian[i02]*_jacobian[i10] -
			 _jacobian[i00]*_jacobian[i12]) / det;
    _jacobianInv[i20] = (_jacobian[i10]*_jacobian[i21] -
			 _jacobian[i11]*_jacobian[i20]) / det;
    _jacobianInv[i21] = (_jacobian[i01]*_jacobian[i20] -
			 _jacobian[i00]*_jacobian[i21]) / det;
    _jacobianInv[i22] = (_jacobian[i00]*_jacobian[i11] -
			 _jacobian[i01]*_jacobian[i10]) / det;

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
  
  PetscLogFlops(numQuadPts*(2+36 + 
			    numBasis*spaceDim*cellDim*4));
} // computeGeometry


// End of file 
