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

#include <math.h> // USES fabs()

#include <assert.h> // USES assert()
#include <stdexcept> // USES internal_error

// ----------------------------------------------------------------------
// Constructor
pylith::feassemble::Quadrature2Din3D::Quadrature2Din3D(void)
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
  Quadrature(q)
{ // copy constructor
} // copy constructor

// ----------------------------------------------------------------------
// Compute geometric quantities for a cell.
void
pylith::feassemble::Quadrature2Din3D::computeGeometry(
		       const ALE::Obj<ALE::Mesh::real_section_type>& coordinates,
		       const ALE::Mesh::point_type& cell)
{ // computeGeometry
  assert(2 == _cellDim);
  assert(3 == _spaceDim);
  assert(0 != _basisDeriv);
  assert(0 != _quadPtsRef);
  assert(0 != _quadPts);
  assert(0 != _quadWts);
  assert(0 != _jacobian);
  assert(0 != _jacobianInv);

  _resetGeometry();

  // Get coordinates of cell's vertices
  const ALE::Mesh::topology_type::patch_type patch  = 0;
  const ALE::Mesh::real_section_type::value_type* vertCoords = 
    coordinates->restrict(patch, cell);
  //assert(3 == coordinates.GetFiberDimensionByDepth(patch,
  //*vertices->begin(), 0));

  // Loop over quadrature points
  for (int iQuadPt=0; iQuadPt < _numQuadPts; ++iQuadPt) {
    
    // Compute coordinates of quadrature point in cell
    // x = sum[i=0,n-1] (Ni * xi)
    // y = sum[i=0,n-1] (Ni * yi)
    // z = sum[i=0,n-1] (Ni * zi)
    for (int iVertex=0, iB=iQuadPt*_numCorners;
	 iVertex < _numCorners;
	 ++iVertex) {
      const double basis = _basis[iB+iVertex];
      for (int iDim=0, iQ=iQuadPt*_spaceDim, iV=iVertex*_spaceDim;
	   iDim < _spaceDim;
	   ++iDim)
	_quadPts[iQ+iDim] +=  basis * vertCoords[iV+iDim];
    } // for
    
    // Compute Jacobian at quadrature point
    // J = [dx/dp dy/dp dz/dp]
    //     [dx/dq dy/dq dz/dq]
    // dx/dp = sum[i=0,n-1] (dNi/dp * xi)
    // dy/dp = sum[i=0,n-1] (dNi/dp * yi)
    // dz/dp = sum[i=0,n-1] (dNi/dp * zi)
    // dx/dq = sum[i=0,n-1] (dNi/dq * xi)
    // dy/dq = sum[i=0,n-1] (dNi/dq * yi)
    // dz/dq = sum[i=0,n-1] (dNi/dq * zi)
    for (int iVertex=0; iVertex < _numCorners; ++iVertex)
      for (int iRow=0, 
	     iB=iQuadPt*_numCorners*_cellDim+iVertex*_cellDim;
	   iRow < _cellDim;
	   ++iRow) {
      const double deriv = _basisDeriv[iB+iRow];
      for (int iCol=0, iJ=iQuadPt*_cellDim*_spaceDim + iRow*_spaceDim;
	     iCol < _spaceDim;
	     ++iCol)
	_jacobian[iJ+iCol] += deriv * vertCoords[iVertex*+_spaceDim+iCol];
    } // for

    // Compute determinant of Jacobian at quadrature point
    // |J| = sqrt(J transpose(J))
    const int iJ = iQuadPt*_cellDim*_spaceDim;
    const int i00 = iJ + 0*_spaceDim + 0;
    const int i01 = iJ + 0*_spaceDim + 1;
    const int i02 = iJ + 0*_spaceDim + 2;
    const int i10 = iJ + 1*_spaceDim + 0;
    const int i11 = iJ + 1*_spaceDim + 1;
    const int i12 = iJ + 1*_spaceDim + 2;
    // JJ = J transpose(J)
    const double jj00 = 
      _jacobian[i00]*_jacobian[i00] +
      _jacobian[i01]*_jacobian[i01] +
      _jacobian[i02]*_jacobian[i02];
    const double jj01 =
      _jacobian[i00]*_jacobian[i10] +
      _jacobian[i01]*_jacobian[i11] +
      _jacobian[i02]*_jacobian[i12];
    const double jj10 = jj01;
    const double jj11 = 
      _jacobian[i10]*_jacobian[i10] +
      _jacobian[i11]*_jacobian[i11] +
      _jacobian[i12]*_jacobian[i12];
    const double det = sqrt(jj00*jj11 - jj01*jj10);
    _checkJacobianDet(det);
    _jacobianDet[iQuadPt] = det;
    
    // Compute inverse of Jacobian at quadrature point
    const double d01 = 
      _jacobian[i00]*_jacobian[i11] - _jacobian[i10]*_jacobian[i01];
    const double d12 = 
      _jacobian[i01]*_jacobian[i12] - _jacobian[i11]*_jacobian[i02];
    const double d02 = 
      _jacobian[i00]*_jacobian[i12] - _jacobian[i10]*_jacobian[i02];
    if (fabs(d01) > _minJacobian) {
      // Jinv00 = 1/d01 * J11
      // Jinv01 = 1/d01 * -J01
      // Jinv10 = 1/d01 * -J10
      // Jinv11 = 1/d01 * J00
      _jacobianInv[iJ+0] = _jacobian[i11] / d01; // Jinv00
      _jacobianInv[iJ+1] = -_jacobian[i01] / d01; // Jinv01
      _jacobianInv[iJ+2] = -_jacobian[i10] / d01; // Jinv10
      _jacobianInv[iJ+3] = _jacobian[i00] / d01; // Jinv11
      if (fabs(d12) > _minJacobian) {
	// Jinv20 = 1/d12 -J11
	// Jinv21 = 1/d12 J01
	_jacobianInv[iJ+4] = -_jacobian[i11] / d12; // Jinv20
	_jacobianInv[iJ+5] = _jacobian[i01] / d12; // Jinv21

      } else if (fabs(d02) > _minJacobian) {
	// Jinv20 = 1/d02 -J10
	// Jinv21 = 1/d02 J00
	_jacobianInv[iJ+4] = -_jacobian[i10] / d02; // Jinv20
	_jacobianInv[iJ+5] = _jacobian[i00] / d02; // Jinv21
      } else {
	_jacobianInv[iJ+4] = 0.0; // Jinv20
	_jacobianInv[iJ+5] = 0.0; // Jinv21
      } // if/else
    } else if (fabs(d02) > _minJacobian) {
      // Jinv00 = 1/d02 * J12
      // Jinv01 = 1/d02 * -J02
      // Jinv20 = 1/d02 * -J10
      // Jinv21 = 1/d02 * J00
      _jacobianInv[iJ+0] = _jacobian[i12] / d02; // Jinv00
      _jacobianInv[iJ+1] = -_jacobian[i02] / d02; // Jinv01
      _jacobianInv[iJ+4] = -_jacobian[i10] / d02; // Jinv20
      _jacobianInv[iJ+5] = _jacobian[i00] / d02; // Jinv21
      if (fabs(d12) > _minJacobian) {
	// Jinv10 = 1/d12 J12
	// Jinv11 = 1/d12 -J02
	_jacobianInv[iJ+2] = -_jacobian[i12] / d12; // Jinv10
	_jacobianInv[iJ+3] = _jacobian[i02] / d12; // Jinv11
      } else {
	_jacobianInv[iJ+2] = 0.0; // Jinv10
	_jacobianInv[iJ+3] = 0.0; // Jinv11
      } // if/else
    } else if (fabs(d12) > _minJacobian) {
      _jacobianInv[iJ+0] = 0.0; // Jinv00
      _jacobianInv[iJ+1] = 0.0; // Jinv01
      // Jinv10 = 1/d12 J12
      // Jinv11 = 1/d12 -J02
      // Jinv20 = 1/d12 -J11
      // Jinv21 = 1/d12 J01
      _jacobianInv[iJ+2] = _jacobian[i12] / d12; // Jinv10
      _jacobianInv[iJ+3] = -_jacobian[i02] / d12; // Jin11
      _jacobianInv[iJ+4] = -_jacobian[i11] / d12; // Jinv20
      _jacobianInv[iJ+5] = _jacobian[i01] / d12; // Jinv21
    } else
      throw std::runtime_error("Could not invert Jacobian.");
  } // for
} // computeGeometry

// End of file 
