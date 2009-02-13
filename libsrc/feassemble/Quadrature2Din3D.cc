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

#include "Quadrature2Din3D.hh" // implementation of class methods

#include "QuadratureRefCell.hh" // USES QuadratureRefCell
#include "CellGeometry.hh" // USES CellGeometry

#include "petsc.h" // USES PetscLogFlops

#include <cmath> // USES fabs()
#include <cassert> // USES assert()
#include <stdexcept> // USES std::runtime_error()

#define ISOPARAMETRIC

// ----------------------------------------------------------------------
// Constructor
pylith::feassemble::Quadrature2Din3D::Quadrature2Din3D(const QuadratureRefCell& q) :
  QuadratureEngine(q)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor
pylith::feassemble::Quadrature2Din3D::~Quadrature2Din3D(void)
{ // destructor
} // destructor
  
// ----------------------------------------------------------------------
// Copy constructor.
pylith::feassemble::Quadrature2Din3D::Quadrature2Din3D(const Quadrature2Din3D& q) :
  QuadratureEngine(q)
{ // copy constructor
} // copy constructor

// ----------------------------------------------------------------------
// Compute geometric quantities for a cell at quadrature points.
void
pylith::feassemble::Quadrature2Din3D::computeGeometry(const double* vertCoords,
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
  const double minJacobian = _quadRefCell.minJacobian();

  assert(2 == cellDim);
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
    geometry.coordsRefToGlobal(&_quadPts[iQuadPt*spaceDim],
			       &quadPtsRef[iQuadPt*cellDim],
			       vertCoords, spaceDim);
#endif

#if defined(ISOPARAMETRIC)
    // Compute Jacobian at quadrature point
    // J = [dx/dp dx/dq]
    //     [dy/dp dy/dq]
    //     [dz/dp dz/dq]
    // dx/dp = sum[i=0,n-1] (dNi/dp * xi)
    // dx/dq = sum[i=0,n-1] (dNi/dq * xi)
    // dy/dp = sum[i=0,n-1] (dNi/dp * yi)
    // dy/dq = sum[i=0,n-1] (dNi/dq * yi)
    // dz/dp = sum[i=0,n-1] (dNi/dp * zi)
    // dz/dq = sum[i=0,n-1] (dNi/dq * zi)
    for (int iBasis=0; iBasis < numBasis; ++iBasis)
      for (int iCol=0; iCol < cellDim; ++iCol) {
	const double deriv = 
	  basisDerivRef[iQuadPt*numBasis*cellDim+iBasis*cellDim+iCol];
	for (int iRow=0; iRow < spaceDim; ++iRow)
	  _jacobian[iQuadPt*cellDim*spaceDim+iRow*cellDim+iCol] +=
	    deriv * vertCoords[iBasis*+spaceDim+iRow];
      } // for
    
    // Compute determinant of Jacobian at quadrature point
    // |J| = sqrt(transpose(J) J)
    const int iJ = iQuadPt*cellDim*spaceDim;
    const int i00 = iJ + 0*cellDim + 0;
    const int i01 = iJ + 0*cellDim + 1;
    const int i10 = iJ + 1*cellDim + 0;
    const int i11 = iJ + 1*cellDim + 1;
    const int i20 = iJ + 2*cellDim + 0;
    const int i21 = iJ + 2*cellDim + 1;
    // JJ = transpose(J) J 
    const double jj00 = 
      _jacobian[i00]*_jacobian[i00] +
      _jacobian[i10]*_jacobian[i10] +
      _jacobian[i20]*_jacobian[i20];
    const double jj10 =
      _jacobian[i00]*_jacobian[i01] +
      _jacobian[i10]*_jacobian[i11] +
      _jacobian[i20]*_jacobian[i21];
    const double jj01 = jj10;
    const double jj11 = 
      _jacobian[i01]*_jacobian[i01] +
      _jacobian[i11]*_jacobian[i11] +
      _jacobian[i21]*_jacobian[i21];
    const double det = sqrt(jj00*jj11 - jj01*jj10);
    _checkJacobianDet(det, cell);
    _jacobianDet[iQuadPt] = det;
#else
    // Compute Jacobian and determinant of Jacobian at quadrature point
    geometry.jacobian(&_jacobian[iQuadPt*cellDim*spaceDim],
		      &_jacobianDet[iQuadPt],
		      vertCoords, &quadPtsRef[iQuadPt*cellDim], spaceDim);
    _checkJacobianDet(_jacobianDet[iQuadPt], cell);

    const int iJ = iQuadPt*cellDim*spaceDim;
    const int i00 = iJ + 0*cellDim + 0;
    const int i01 = iJ + 0*cellDim + 1;
    const int i10 = iJ + 1*cellDim + 0;
    const int i11 = iJ + 1*cellDim + 1;
    const int i20 = iJ + 2*cellDim + 0;
    const int i21 = iJ + 2*cellDim + 1;
#endif
    
    // Compute inverse of Jacobian at quadrature point
    const double d01 = 
      _jacobian[i00]*_jacobian[i11] - 
      _jacobian[i10]*_jacobian[i01];
    const double d12 = 
      _jacobian[i10]*_jacobian[i21] - 
      _jacobian[i20]*_jacobian[i11];
    const double d02 = 
      _jacobian[i00]*_jacobian[i21] - 
      _jacobian[i20]*_jacobian[i01];
    if (fabs(d01) > minJacobian) {
      // Jinv00 = 1/d01 * J11
      // Jinv01 = 1/d01 * -J01
      // Jinv10 = 1/d01 * -J10
      // Jinv11 = 1/d01 * J00
      _jacobianInv[iJ+0] =  _jacobian[i11] / d01; // Jinv00
      _jacobianInv[iJ+1] = -_jacobian[i01] / d01; // Jinv01
      _jacobianInv[iJ+3] = -_jacobian[i10] / d01; // Jinv10
      _jacobianInv[iJ+4] =  _jacobian[i00] / d01; // Jinv11
      if (fabs(d12) > minJacobian) {
	// Jinv02 = 1/d12 -J11
	// Jinv12 = 1/d12 J10
	_jacobianInv[iJ+2] = -_jacobian[i11] / d12; // Jinv02
	_jacobianInv[iJ+5] =  _jacobian[i10] / d12; // Jinv12
	
      } else if (fabs(d02) > minJacobian) {
	// Jinv02 = 1/d02 -J01
	// Jinv12 = 1/d02 J00
	_jacobianInv[iJ+2] = -_jacobian[i01] / d02; // Jinv02
	_jacobianInv[iJ+5] =  _jacobian[i00] / d02; // Jinv12
      } else {
	_jacobianInv[iJ+2] = 0.0; // Jinv02
	_jacobianInv[iJ+5] = 0.0; // Jinv12
      } // if/else
    } else if (fabs(d02) > minJacobian) {
      // Jinv00 = 1/d02 * J21
      // Jinv02 = 1/d02 * -J01
      // Jinv10 = 1/d02 * -J20
      // Jinv12 = 1/d02 * J00
      _jacobianInv[iJ+0] =  _jacobian[i21] / d02; // Jinv00
      _jacobianInv[iJ+2] = -_jacobian[i01] / d02; // Jinv02
      _jacobianInv[iJ+3] = -_jacobian[i20] / d02; // Jinv10
      _jacobianInv[iJ+5] =  _jacobian[i00] / d02; // Jinv12
      if (fabs(d12) > minJacobian) {
	// Jinv01 = 1/d12 J21
	// Jinv11 = 1/d12 -J20
	_jacobianInv[iJ+1] = -_jacobian[i21] / d12; // Jinv01
	_jacobianInv[iJ+4] =  _jacobian[i20] / d12; // Jinv11
      } else {
	_jacobianInv[iJ+1] = 0.0; // Jinv01
	_jacobianInv[iJ+4] = 0.0; // Jinv11
      } // if/else
    } else if (fabs(d12) > minJacobian) {
      _jacobianInv[iJ+0] = 0.0; // Jinv00
      _jacobianInv[iJ+3] = 0.0; // Jinv10
      // Jinv01 = 1/d12 J21
      // Jinv02 = 1/d12 -J11
      // Jinv11 = 1/d12 -J20
      // Jinv12 = 1/d12 J10
      _jacobianInv[iJ+1] =  _jacobian[i21] / d12; // Jinv01
      _jacobianInv[iJ+2] = -_jacobian[i11] / d12; // Jinv02
      _jacobianInv[iJ+4] = -_jacobian[i20] / d12; // Jinv11
      _jacobianInv[iJ+5] =  _jacobian[i10] / d12; // Jinv12
    } else
      throw std::runtime_error("Could not invert Jacobian.");

    // Compute derivatives of basis functions with respect to global
    // coordinates
    // dNi/dx = dNi/dp dp/dx + dNi/dq dq/dx + dNi/dr dr/dx
    for (int iBasis=0; iBasis < numBasis; ++iBasis)
      for (int iDim=0; iDim < spaceDim; ++iDim)
	for (int jDim=0; jDim < cellDim; ++jDim)
	  _basisDeriv[iQuadPt*numBasis*spaceDim+iBasis*spaceDim+iDim] +=
	    basisDerivRef[iQuadPt*numBasis*cellDim + iBasis*cellDim+jDim] *
	    _jacobianInv[iQuadPt*cellDim*spaceDim+jDim*spaceDim+iDim];
  } // for
  
  PetscLogFlops(numQuadPts*(15 +
			    numBasis*spaceDim*2 +
			    numBasis*spaceDim*cellDim*2));
} // computeGeometry


// End of file 
