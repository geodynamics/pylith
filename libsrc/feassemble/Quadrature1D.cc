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
pylith::feassemble::Quadrature1D::computeGeometry(const double* vertCoords,
						  const int coordDim,
						  const int cell)
{ // computeGeometry
  assert(0 != vertCoords);
  assert(1 == coordDim);

  const int cellDim = _quadRefCell.cellDim();
  const int spaceDim = _quadRefCell.spaceDim();
  const int numQuadPts = _quadRefCell.numQuadPts();
  const int numBasis = _quadRefCell.numBasis();

  const double_array& basis = _quadRefCell.basis();
  const double_array& quadPtsRef = _quadRefCell.quadPtsRef();
  const double_array& basisDerivRef = _quadRefCell.basisDerivRef();
  const CellGeometry& geometry = _quadRefCell.refGeometry();

  assert(1 == cellDim);
  assert(1 == spaceDim);
  zero();
  
  // Loop over quadrature points
  for (int iQuadPt=0; iQuadPt < numQuadPts; ++iQuadPt) {

    // Compute coordinates of quadrature point in cell
#if defined(ISOPARAMETRIC)
    // x = sum[i=0,n-1] (Ni * xi)
    for (int iBasis=0; iBasis < numBasis; ++iBasis)
      _quadPts[iQuadPt] += 
	basis[iQuadPt*numBasis+iBasis]*vertCoords[iBasis];
#else
    geometry.coordsRefToGlobal(&_quadPts[iQuadPt], &quadPtsRef[iQuadPt],
			       vertCoords, spaceDim);
#endif

#if defined(ISOPARAMETRIC)
    // Compute Jacobian at quadrature point
    // J = dx/dp = sum[i=0,n-1] (dNi/dp * xi)
    for (int iBasis=0; iBasis < numBasis; ++iBasis)
      _jacobian[iQuadPt] += 
	basisDerivRef[iQuadPt*numBasis+iBasis] * vertCoords[iBasis];

    // Compute determinant of Jacobian at quadrature point
    // |J| = j00
    const double det = _jacobian[iQuadPt];
    _checkJacobianDet(det, cell);
    _jacobianDet[iQuadPt] = _jacobian[iQuadPt];
#else
    // Compute Jacobian and determinant of Jacobian at quadrature point
    assert(0 != _geometry);
    geometry->jacobian(&_jacobian[iQuadPt], &_jacobianDet[iQuadPt],
		       vertCoords, &quadPtsRef[iQuadPt], spaceDim);
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
      _basisDeriv[iQuadPt*numBasis+iBasis] +=
	  basisDerivRef[iQuadPt*numBasis+iBasis] *
	  _jacobianInv[iQuadPt];
  } // for

  PetscLogFlops(numQuadPts * (1 + numBasis * 4));
} // computeGeometry


// End of file 
