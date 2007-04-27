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

#include <assert.h> // USES assert()

// ----------------------------------------------------------------------
// Constructor
pylith::feassemble::Quadrature3D::Quadrature3D(void)
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
  Quadrature(q)
{ // copy constructor
} // copy constructor

// ----------------------------------------------------------------------
// Compute geometric quantities for a cell at vertices.
void
pylith::feassemble::Quadrature3D::computeGeometryVert(
		       const ALE::Obj<Mesh>& mesh,
		       const ALE::Obj<real_section_type>& coordinates,
		       const Mesh::point_type& cell)
{ // computeGeometryVert
  assert(3 == _cellDim);
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
    //     [dx/dr dy/dr dz/dr]
    // dx/dp = sum[i=0,n-1] (dNi/dp * xi)
    // dy/dp = sum[i=0,n-1] (dNi/dp * yi)
    // dz/dp = sum[i=0,n-1] (dNi/dp * zi)
    // dx/dq = sum[i=0,n-1] (dNi/dq * xi)
    // dy/dq = sum[i=0,n-1] (dNi/dq * yi)
    // dz/dq = sum[i=0,n-1] (dNi/dq * zi)
    // dx/dr = sum[i=0,n-1] (dNi/dr * xi)
    // dy/dr = sum[i=0,n-1] (dNi/dr * yi)
    // dz/dr = sum[i=0,n-1] (dNi/dr * zi)
    for (int iBasis=0; iBasis < _numBasis; ++iBasis)
      for (int iRow=0; iRow < _cellDim; ++iRow) {
	const double deriv = 
	  _basisDerivVert[iVertex*_numBasis*_spaceDim+iBasis*_cellDim+iRow];
	for (int iCol=0; iCol < _spaceDim; ++iCol)
	  _jacobianVert[iVertex*_cellDim*_spaceDim+iRow*_spaceDim+iCol] += 
	    deriv * vertCoords[iBasis*_spaceDim+iCol];
      } // for

    // Compute determinant of Jacobian at vertex
    // |J| = j00*(j11*j22-j12*j21) +
    //      -j01*(j10*j22-j12*j20) +
    //       j02*(j10*j21-j11*j20)
    const int iJ = iVertex*_cellDim*_spaceDim;
    const int i00 = iJ + 0*_spaceDim + 0;
    const int i01 = iJ + 0*_spaceDim + 1;
    const int i02 = iJ + 0*_spaceDim + 2;
    const int i10 = iJ + 1*_spaceDim + 0;
    const int i11 = iJ + 1*_spaceDim + 1;
    const int i12 = iJ + 1*_spaceDim + 2;
    const int i20 = iJ + 2*_spaceDim + 0;
    const int i21 = iJ + 2*_spaceDim + 1;
    const int i22 = iJ + 2*_spaceDim + 2;
    const double det = 
      _jacobianVert[i00]*(_jacobianVert[i11]*_jacobianVert[i22] -
		      _jacobianVert[i12]*_jacobianVert[i21]) -
      _jacobianVert[i01]*(_jacobianVert[i10]*_jacobianVert[i22] -
		      _jacobianVert[i12]*_jacobianVert[i20]) +
      _jacobianVert[i02]*(_jacobianVert[i10]*_jacobianVert[i21] -
		      _jacobianVert[i11]*_jacobianVert[i20]);
    _checkJacobianDet(det);
    _jacobianDetVert[iVertex] = det;
  } // for
    
} // computeGeometryVert

// ----------------------------------------------------------------------
// Compute geometric quantities for a cell at quadrature points.
void
pylith::feassemble::Quadrature3D::computeGeometryQuad(
		       const ALE::Obj<Mesh>& mesh,
		       const ALE::Obj<real_section_type>& coordinates,
		       const Mesh::point_type& cell)
{ // computeGeometryQuad
  assert(3 == _cellDim);
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
    //     [dx/dr dy/dr dz/dr]
    // dx/dp = sum[i=0,n-1] (dNi/dp * xi)
    // dy/dp = sum[i=0,n-1] (dNi/dp * yi)
    // dz/dp = sum[i=0,n-1] (dNi/dp * zi)
    // dx/dq = sum[i=0,n-1] (dNi/dq * xi)
    // dy/dq = sum[i=0,n-1] (dNi/dq * yi)
    // dz/dq = sum[i=0,n-1] (dNi/dq * zi)
    // dx/dr = sum[i=0,n-1] (dNi/dr * xi)
    // dy/dr = sum[i=0,n-1] (dNi/dr * yi)
    // dz/dr = sum[i=0,n-1] (dNi/dr * zi)
    for (int iBasis=0; iBasis < _numBasis; ++iBasis)
      for (int iRow=0; iRow < _cellDim; ++iRow) {
	const double deriv = 
	  _basisDerivQuad[iQuadPt*_numBasis*_spaceDim+iBasis*_cellDim+iRow];
	for (int iCol=0; iCol < _spaceDim; ++iCol)
	  _jacobianQuad[iQuadPt*_cellDim*_spaceDim+iRow*_spaceDim+iCol] += 
	    deriv * vertCoords[iBasis*_spaceDim+iCol];
      } // for

    // Compute determinant of Jacobian at quadrature point
    // |J| = j00*(j11*j22-j12*j21) +
    //      -j01*(j10*j22-j12*j20) +
    //       j02*(j10*j21-j11*j20)
    const int iJ = iQuadPt*_cellDim*_spaceDim;
    const int i00 = iJ + 0*_spaceDim + 0;
    const int i01 = iJ + 0*_spaceDim + 1;
    const int i02 = iJ + 0*_spaceDim + 2;
    const int i10 = iJ + 1*_spaceDim + 0;
    const int i11 = iJ + 1*_spaceDim + 1;
    const int i12 = iJ + 1*_spaceDim + 2;
    const int i20 = iJ + 2*_spaceDim + 0;
    const int i21 = iJ + 2*_spaceDim + 1;
    const int i22 = iJ + 2*_spaceDim + 2;
    const double det = 
      _jacobianQuad[i00]*(_jacobianQuad[i11]*_jacobianQuad[i22] -
		      _jacobianQuad[i12]*_jacobianQuad[i21]) -
      _jacobianQuad[i01]*(_jacobianQuad[i10]*_jacobianQuad[i22] -
		      _jacobianQuad[i12]*_jacobianQuad[i20]) +
      _jacobianQuad[i02]*(_jacobianQuad[i10]*_jacobianQuad[i21] -
		      _jacobianQuad[i11]*_jacobianQuad[i20]);
    _checkJacobianDet(det);
    _jacobianDetQuad[iQuadPt] = det;
    
    // Compute inverse of Jacobian at quadrature point
    _jacobianInvQuad[i00] = (_jacobianQuad[i11]*_jacobianQuad[i22] -
			 _jacobianQuad[i12]*_jacobianQuad[i21]) / det;
    _jacobianInvQuad[i01] = (_jacobianQuad[i02]*_jacobianQuad[i21] -
			 _jacobianQuad[i01]*_jacobianQuad[i22]) / det;
    _jacobianInvQuad[i02] = (_jacobianQuad[i01]*_jacobianQuad[i12] -
			 _jacobianQuad[i02]*_jacobianQuad[i11]) / det;
    _jacobianInvQuad[i10] = (_jacobianQuad[i12]*_jacobianQuad[i20] -
			 _jacobianQuad[i10]*_jacobianQuad[i22]) / det;
    _jacobianInvQuad[i11] = (_jacobianQuad[i00]*_jacobianQuad[i22] -
			 _jacobianQuad[i02]*_jacobianQuad[i20]) / det;
    _jacobianInvQuad[i12] = (_jacobianQuad[i02]*_jacobianQuad[i10] -
			 _jacobianQuad[i00]*_jacobianQuad[i12]) / det;
    _jacobianInvQuad[i20] = (_jacobianQuad[i10]*_jacobianQuad[i21] -
			 _jacobianQuad[i11]*_jacobianQuad[i20]) / det;
    _jacobianInvQuad[i21] = (_jacobianQuad[i01]*_jacobianQuad[i20] -
			 _jacobianQuad[i00]*_jacobianQuad[i21]) / det;
    _jacobianInvQuad[i22] = (_jacobianQuad[i00]*_jacobianQuad[i11] -
			 _jacobianQuad[i01]*_jacobianQuad[i10]) / det;
  } // for
} // computeGeometryQuad


// End of file 
