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
// Compute geometric quantities for a cell at vertices.
void
pylith::feassemble::Quadrature2Din3D::computeGeometryVert(
		       const ALE::Obj<Mesh>& mesh,
		       const ALE::Obj<real_section_type>& coordinates,
		       const Mesh::point_type& cell)
{ // computeGeometryVert
  assert(2 == _cellDim);
  assert(3 == _spaceDim);

  _resetGeometry();

  // Get coordinates of cell's vertices
  const real_section_type::value_type* vertCoords = 
    mesh->restrict(coordinates, cell);
  assert(3 == coordinates->getFiberDimension(*mesh->depthStratum(0)->begin()));

  // Loop over vertices
  const int numVertices = _numBasis;
  for (int iVertex=0; iVertex < numVertices; ++iVertex) {
    // Compute Jacobian at vertex
    // J = [dx/dp dy/dp dz/dp]
    //     [dx/dq dy/dq dz/dq]
    // dx/dp = sum[i=0,n-1] (dNi/dp * xi)
    // dy/dp = sum[i=0,n-1] (dNi/dp * yi)
    // dz/dp = sum[i=0,n-1] (dNi/dp * zi)
    // dx/dq = sum[i=0,n-1] (dNi/dq * xi)
    // dy/dq = sum[i=0,n-1] (dNi/dq * yi)
    // dz/dq = sum[i=0,n-1] (dNi/dq * zi)
    for (int iBasis=0; iBasis < _numBasis; ++iBasis)
      for (int iRow=0; iRow < _cellDim; ++iRow) {
	const double deriv = 
	  _basisDerivVert[iVertex*_numBasis*_cellDim+iBasis*_cellDim+iRow];
	for (int iCol=0; iCol < _spaceDim; ++iCol)
	  _jacobianVert[iVertex*_cellDim*_spaceDim+iRow*_spaceDim+iCol] +=
	    deriv * vertCoords[iBasis*+_spaceDim+iCol];
      } // for
    
    // Compute determinant of Jacobian at quadrature point
    // |J| = sqrt(J transpose(J))
    const int iJ = iVertex*_cellDim*_spaceDim;
    const int i00 = iJ + 0*_spaceDim + 0;
    const int i01 = iJ + 0*_spaceDim + 1;
    const int i02 = iJ + 0*_spaceDim + 2;
    const int i10 = iJ + 1*_spaceDim + 0;
    const int i11 = iJ + 1*_spaceDim + 1;
    const int i12 = iJ + 1*_spaceDim + 2;
    // JJ = J transpose(J)
    const double jj00 = 
      _jacobianVert[i00]*_jacobianVert[i00] +
      _jacobianVert[i01]*_jacobianVert[i01] +
      _jacobianVert[i02]*_jacobianVert[i02];
    const double jj01 =
      _jacobianVert[i00]*_jacobianVert[i10] +
      _jacobianVert[i01]*_jacobianVert[i11] +
      _jacobianVert[i02]*_jacobianVert[i12];
    const double jj10 = jj01;
    const double jj11 = 
      _jacobianVert[i10]*_jacobianVert[i10] +
      _jacobianVert[i11]*_jacobianVert[i11] +
      _jacobianVert[i12]*_jacobianVert[i12];
    const double det = sqrt(jj00*jj11 - jj01*jj10);
    _checkJacobianDet(det);
    _jacobianDetVert[iVertex] = det;
  } // for
    
} // computeGeometryVert

// ----------------------------------------------------------------------
// Compute geometric quantities for a cell at quadrature points.
void
pylith::feassemble::Quadrature2Din3D::computeGeometryQuad(
		       const ALE::Obj<Mesh>& mesh,
		       const ALE::Obj<real_section_type>& coordinates,
		       const Mesh::point_type& cell)
{ // computeGeometryQuad
  assert(2 == _cellDim);
  assert(3 == _spaceDim);

  _resetGeometry();

  // Get coordinates of cell's vertices
  const real_section_type::value_type* vertCoords = 
    mesh->restrict(coordinates, cell);
  assert(3 == coordinates->getFiberDimension(*mesh->depthStratum(0)->begin()));

  // Loop over quadrature points
  for (int iQuadPt=0; iQuadPt < _numQuadPts; ++iQuadPt) {
    
    // Compute coordinates of quadrature point in cell
    // x = sum[i=0,n-1] (Ni * xi)
    // y = sum[i=0,n-1] (Ni * yi)
    // z = sum[i=0,n-1] (Ni * zi)
    for (int iBasis=0; iBasis < _numBasis; ++iBasis) {
      const double basis = _basisQuad[iQuadPt*_numBasis+iBasis];
      for (int iDim=0; iDim < _spaceDim; ++iDim)
	_quadPts[iQuadPt*_spaceDim+iDim] += 
	  basis * vertCoords[iBasis*_spaceDim+iDim];
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
    for (int iBasis=0; iBasis < _numBasis; ++iBasis)
      for (int iRow=0; iRow < _cellDim; ++iRow) {
	const double deriv = 
	  _basisDerivQuad[iQuadPt*_numBasis*_cellDim+iBasis*_cellDim+iRow];
	for (int iCol=0; iCol < _spaceDim; ++iCol)
	  _jacobianQuad[iQuadPt*_cellDim*_spaceDim+iRow*_spaceDim+iCol] +=
	    deriv * vertCoords[iBasis*+_spaceDim+iCol];
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
      _jacobianQuad[i00]*_jacobianQuad[i00] +
      _jacobianQuad[i01]*_jacobianQuad[i01] +
      _jacobianQuad[i02]*_jacobianQuad[i02];
    const double jj01 =
      _jacobianQuad[i00]*_jacobianQuad[i10] +
      _jacobianQuad[i01]*_jacobianQuad[i11] +
      _jacobianQuad[i02]*_jacobianQuad[i12];
    const double jj10 = jj01;
    const double jj11 = 
      _jacobianQuad[i10]*_jacobianQuad[i10] +
      _jacobianQuad[i11]*_jacobianQuad[i11] +
      _jacobianQuad[i12]*_jacobianQuad[i12];
    const double det = sqrt(jj00*jj11 - jj01*jj10);
    _checkJacobianDet(det);
    _jacobianDetQuad[iQuadPt] = det;
    
    // Compute inverse of Jacobian at quadrature point
    const double d01 = 
      _jacobianQuad[i00]*_jacobianQuad[i11] - 
      _jacobianQuad[i10]*_jacobianQuad[i01];
    const double d12 = 
      _jacobianQuad[i01]*_jacobianQuad[i12] - 
      _jacobianQuad[i11]*_jacobianQuad[i02];
    const double d02 = 
      _jacobianQuad[i00]*_jacobianQuad[i12] - 
      _jacobianQuad[i10]*_jacobianQuad[i02];
    if (fabs(d01) > _minJacobian) {
      // Jinv00 = 1/d01 * J11
      // Jinv01 = 1/d01 * -J01
      // Jinv10 = 1/d01 * -J10
      // Jinv11 = 1/d01 * J00
      _jacobianInvQuad[iJ+0] =  _jacobianQuad[i11] / d01; // Jinv00
      _jacobianInvQuad[iJ+1] = -_jacobianQuad[i01] / d01; // Jinv01
      _jacobianInvQuad[iJ+2] = -_jacobianQuad[i10] / d01; // Jinv10
      _jacobianInvQuad[iJ+3] =  _jacobianQuad[i00] / d01; // Jinv11
      if (fabs(d12) > _minJacobian) {
	// Jinv20 = 1/d12 -J11
	// Jinv21 = 1/d12 J01
	_jacobianInvQuad[iJ+4] = -_jacobianQuad[i11] / d12; // Jinv20
	_jacobianInvQuad[iJ+5] =  _jacobianQuad[i01] / d12; // Jinv21
	
      } else if (fabs(d02) > _minJacobian) {
	// Jinv20 = 1/d02 -J10
	// Jinv21 = 1/d02 J00
	_jacobianInvQuad[iJ+4] = -_jacobianQuad[i10] / d02; // Jinv20
	_jacobianInvQuad[iJ+5] =  _jacobianQuad[i00] / d02; // Jinv21
      } else {
	_jacobianInvQuad[iJ+4] = 0.0; // Jinv20
	_jacobianInvQuad[iJ+5] = 0.0; // Jinv21
      } // if/else
    } else if (fabs(d02) > _minJacobian) {
      // Jinv00 = 1/d02 * J12
      // Jinv01 = 1/d02 * -J02
      // Jinv20 = 1/d02 * -J10
      // Jinv21 = 1/d02 * J00
      _jacobianInvQuad[iJ+0] =  _jacobianQuad[i12] / d02; // Jinv00
      _jacobianInvQuad[iJ+1] = -_jacobianQuad[i02] / d02; // Jinv01
      _jacobianInvQuad[iJ+4] = -_jacobianQuad[i10] / d02; // Jinv20
      _jacobianInvQuad[iJ+5] =  _jacobianQuad[i00] / d02; // Jinv21
      if (fabs(d12) > _minJacobian) {
	// Jinv10 = 1/d12 J12
	// Jinv11 = 1/d12 -J02
	_jacobianInvQuad[iJ+2] = -_jacobianQuad[i12] / d12; // Jinv10
	_jacobianInvQuad[iJ+3] =  _jacobianQuad[i02] / d12; // Jinv11
      } else {
	_jacobianInvQuad[iJ+2] = 0.0; // Jinv10
	_jacobianInvQuad[iJ+3] = 0.0; // Jinv11
      } // if/else
    } else if (fabs(d12) > _minJacobian) {
      _jacobianInvQuad[iJ+0] = 0.0; // Jinv00
      _jacobianInvQuad[iJ+1] = 0.0; // Jinv01
      // Jinv10 = 1/d12 J12
      // Jinv11 = 1/d12 -J02
      // Jinv20 = 1/d12 -J11
      // Jinv21 = 1/d12 J01
      _jacobianInvQuad[iJ+2] =  _jacobianQuad[i12] / d12; // Jinv10
      _jacobianInvQuad[iJ+3] = -_jacobianQuad[i02] / d12; // Jin11
      _jacobianInvQuad[iJ+4] = -_jacobianQuad[i11] / d12; // Jinv20
      _jacobianInvQuad[iJ+5] =  _jacobianQuad[i01] / d12; // Jinv21
    } else
      throw std::runtime_error("Could not invert Jacobian.");
  } // for
} // computeGeometryQuad


// End of file 
