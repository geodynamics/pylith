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
pylith::feassemble::Quadrature2Din3D::computeGeometry(const PylithScalar* coordinatesCell,
						      const int coordinatesSize,
						      const int cell)
{ // computeGeometry
  assert(coordinatesCell);

  const int cellDim = 2;
  const int spaceDim = 3;

  const int numQuadPts = _quadRefCell.numQuadPts();
  const int numBasis = _quadRefCell.numBasis();

  const scalar_array& basis = _quadRefCell.basis();
  const scalar_array& basisDerivRef = _quadRefCell.basisDerivRef();
  const PylithScalar minJacobian = _quadRefCell.minJacobian();

  assert(_quadRefCell.cellDim() == cellDim);
  assert(_quadRefCell.spaceDim() == spaceDim);
  assert(numBasis*spaceDim == coordinatesSize);

  zero();
  
  // Loop over quadrature points
  for (int iQuadPt=0; iQuadPt < numQuadPts; ++iQuadPt) {
    
    // Compute coordinates of quadrature point in cell
#if defined(ISOPARAMETRIC)
    // x = sum[i=0,n-1] (Ni * xi)
    // y = sum[i=0,n-1] (Ni * yi)
    // z = sum[i=0,n-1] (Ni * zi)
    for (int iBasis=0; iBasis < numBasis; ++iBasis) {
      const PylithScalar valueBasis = basis[iQuadPt*numBasis+iBasis];
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
	const PylithScalar deriv = 
	  basisDerivRef[iQuadPt*numBasis*cellDim+iBasis*cellDim+iCol];
	for (int iRow=0; iRow < spaceDim; ++iRow)
	  _jacobian[iQuadPt*cellDim*spaceDim+iRow*cellDim+iCol] +=
	    deriv * coordinatesCell[iBasis*+spaceDim+iRow];
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
    const PylithScalar jj00 = 
      _jacobian[i00]*_jacobian[i00] +
      _jacobian[i10]*_jacobian[i10] +
      _jacobian[i20]*_jacobian[i20];
    const PylithScalar jj10 =
      _jacobian[i00]*_jacobian[i01] +
      _jacobian[i10]*_jacobian[i11] +
      _jacobian[i20]*_jacobian[i21];
    const PylithScalar jj01 = jj10;
    const PylithScalar jj11 = 
      _jacobian[i01]*_jacobian[i01] +
      _jacobian[i11]*_jacobian[i11] +
      _jacobian[i21]*_jacobian[i21];
    const PylithScalar det = sqrt(jj00*jj11 - jj01*jj10);
    _checkJacobianDet(det, cell);
    _jacobianDet[iQuadPt] = det;
#else
    // Compute Jacobian and determinant of Jacobian at quadrature point
    const int iJ = iQuadPt*cellDim*spaceDim;
    const int i00 = iJ + 0*cellDim + 0;
    const int i01 = iJ + 0*cellDim + 1;
    const int i10 = iJ + 1*cellDim + 0;
    const int i11 = iJ + 1*cellDim + 1;
    const int i20 = iJ + 2*cellDim + 0;
    const int i21 = iJ + 2*cellDim + 1;
    geometry.jacobian(&_jacobian[iQuadPt*cellDim*spaceDim],
		      &_jacobianDet[iQuadPt],
		      coordinatesCell, &quadPtsRef[iQuadPt*cellDim], 
		      spaceDim, 1);
    _checkJacobianDet(_jacobianDet[iQuadPt], cell);
#endif
    
    // Compute inverse of Jacobian at quadrature point
    const PylithScalar d01 = 
      _jacobian[i00]*_jacobian[i11] - 
      _jacobian[i10]*_jacobian[i01];
    const PylithScalar d12 = 
      _jacobian[i10]*_jacobian[i21] - 
      _jacobian[i20]*_jacobian[i11];
    const PylithScalar d02 = 
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
    for (int iBasis=0; iBasis < numBasis; ++iBasis) {
      const int iD = iQuadPt*numBasis*spaceDim + iBasis*spaceDim;
      const int iDR = iQuadPt*numBasis*cellDim + iBasis*cellDim;
      for (int iDim=0; iDim < spaceDim; ++iDim)
        for (int jDim=0; jDim < cellDim; ++jDim)
          _basisDeriv[iD+iDim] += basisDerivRef[iDR+jDim] *
              _jacobianInv[iJ+jDim*spaceDim+iDim];
    } // for
  } // for
  
  PetscLogFlops(numQuadPts*(15 +
			    numBasis*spaceDim*2 +
			    numBasis*spaceDim*cellDim*2));
} // computeGeometry


// End of file 
