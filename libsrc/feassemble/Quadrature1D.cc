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

//#include "Quadrature1D.hh" // implementation of class methods

#include "CellGeometry.hh" // USES CellGeometry

#include "petsc.h" // USES PetscLogFlops

#include <cassert> // USES assert()

#define ISOPARAMETRIC

// ----------------------------------------------------------------------
// Constructor
template<typename mesh_type>
pylith::feassemble::Quadrature1D<mesh_type>::Quadrature1D(void)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor
template<typename mesh_type>
pylith::feassemble::Quadrature1D<mesh_type>::~Quadrature1D(void)
{ // destructor
} // destructor
  
// ----------------------------------------------------------------------
// Copy constructor.
template<typename mesh_type>
pylith::feassemble::Quadrature1D<mesh_type>::Quadrature1D(const Quadrature1D& q) :
  Quadrature<mesh_type>(q)
{ // copy constructor
} // copy constructor

#if 1
// ----------------------------------------------------------------------
// Compute geometric quantities for a cell at quadrature points.
template<typename mesh_type>
void
pylith::feassemble::Quadrature1D<mesh_type>::_computeGeometry(
					      const double* vertCoords,
					      const int coordDim,
					      const int cell)
{ // computeGeometry
  const int cellDim = this->_cellDim;
  const int spaceDim = this->_spaceDim;

  assert(1 == this->_cellDim);
  assert(1 == this->_spaceDim);
  assert(1 == coordDim);
  
  // Loop over quadrature points
  for (int iQuadPt=0; iQuadPt < _numQuadPts; ++iQuadPt) {

    // Compute coordinates of quadrature point in cell
#if defined(ISOPARAMETRIC)
    // x = sum[i=0,n-1] (Ni * xi)
    for (int iBasis=0; iBasis < _numBasis; ++iBasis)
      _quadPts[iQuadPt] += 
	_basis[iQuadPt*_numBasis+iBasis]*vertCoords[iBasis];
#else
    assert(0 != _geometry);
    _geometry->coordsRefToGlobal(&_quadPts[iQuadPt], &_quadPtsRef[iQuadPt],
				 vertCoords, _spaceDim);
#endif

#if defined(ISOPARAMETRIC)
    // Compute Jacobian at quadrature point
    // J = dx/dp = sum[i=0,n-1] (dNi/dp * xi)
    for (int iBasis=0; iBasis < _numBasis; ++iBasis)
      _jacobian[iQuadPt] += 
	_basisDerivRef[iQuadPt*_numBasis+iBasis] * vertCoords[iBasis];

    // Compute determinant of Jacobian at quadrature point
    // |J| = j00
    const double det = _jacobian[iQuadPt];
    _checkJacobianDet(det, cell);
    _jacobianDet[iQuadPt] = _jacobian[iQuadPt];
#else
    // Compute Jacobian and determinant of Jacobian at quadrature point
    assert(0 != _geometry);
    _geometry->jacobian(&_jacobian[iQuadPt], &_jacobianDet[iQuadPt],
			vertCoords, &_quadPtsRef[iQuadPt], _spaceDim);
    _checkJacobianDet(_jacobianDet[iQuadPt], cell);
#endif
    
    // Compute inverse of Jacobian at quadrature point
    // Jinv = 1/j00
    _jacobianInv[iQuadPt] = 1.0 / _jacobianDet[iQuadPt];

    assert(_numQuadPts*_numBasis*_spaceDim == _basisDeriv.size());
    assert(_numQuadPts*_numBasis*_cellDim == _basisDerivRef.size());
    assert(_numQuadPts*_cellDim*_spaceDim == _jacobianInv.size());

    // Compute derivatives of basis functions with respect to global
    // coordinates
    // dNi/dx = dNi/dp dp/dx + dNi/dq dq/dx + dNi/dr dr/dx
    for (int iBasis=0; iBasis < _numBasis; ++iBasis)
      _basisDeriv[iQuadPt*_numBasis+iBasis] +=
	  _basisDerivRef[iQuadPt*_numBasis+iBasis] *
	  _jacobianInv[iQuadPt];
  } // for

  PetscLogFlops(_numQuadPts * (1 + _numBasis * 4));
} // computeGeometry
#endif

// End of file 
