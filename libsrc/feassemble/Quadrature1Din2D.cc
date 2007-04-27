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

#include <assert.h> // USES assert()

// ----------------------------------------------------------------------
// Constructor
pylith::feassemble::Quadrature1Din2D::Quadrature1Din2D(void)
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
  Quadrature(q)
{ // copy constructor
} // copy constructor

// ----------------------------------------------------------------------
// Compute geometric quantities for a cell at vertices.
void
pylith::feassemble::Quadrature1Din2D::computeGeometryVert(
		       const ALE::Obj<Mesh>& mesh,
		       const ALE::Obj<real_section_type>& coordinates,
		       const Mesh::point_type& cell)
{ // computeGeometryVert
  assert(1 == _cellDim);
  assert(2 == _spaceDim);

  _resetGeometry();

  // Get coordinates of cell's vertices
  const real_section_type::value_type* vertCoords = 
    mesh->restrict(coordinates, cell);
  assert(2 == coordinates->getFiberDimension(*mesh->depthStratum(0)->begin()));

  // Loop over vertices
  const int numVertices = _numBasis;
  for (int iVertex=0; iVertex < numVertices; ++iVertex) {
    
    // Compute Jacobian at vertex
    // J = [dx/dp dy/dp]
    // dx/dp = sum[i=0,n-1] (dNi/dp * xi)
    // dy/dp = sum[i=0,n-1] (dNi/dp * yi)
    for (int iBasis=0; iBasis < _numBasis; ++iBasis) {
      const double deriv = _basisDerivVert[iVertex*_numBasis+iBasis];
      for (int iDim=0; iDim < _spaceDim; ++iDim)
	_jacobianVert[iVertex*_spaceDim+iDim] += 
	  deriv * vertCoords[iBasis*_spaceDim+iDim];
    } // for

    // Compute determinant of Jacobian at vertex
    // |J| = sqrt(J transpose(J))
    double det = 0.0;
    for (int iDim=0; iDim < _spaceDim; ++iDim)
      det += _jacobianVert[iVertex*_spaceDim+iDim] * 
	_jacobianVert[iVertex*_spaceDim+iDim];
    det = sqrt(det);
    _checkJacobianDet(det);
    _jacobianDetVert[iVertex] = det;
  } // for

} // computeGeometryVert

// ----------------------------------------------------------------------
// Compute geometric quantities for a cell at quadrature points.
void
pylith::feassemble::Quadrature1Din2D::computeGeometryQuad(
		       const ALE::Obj<Mesh>& mesh,
		       const ALE::Obj<real_section_type>& coordinates,
		       const Mesh::point_type& cell)
{ // computeGeometryQuad
  assert(1 == _cellDim);
  assert(2 == _spaceDim);

  _resetGeometry();

  // Get coordinates of cell's vertices
  const real_section_type::value_type* vertCoords = 
    mesh->restrict(coordinates, cell);
  assert(2 == coordinates->getFiberDimension(*mesh->depthStratum(0)->begin()));

  // Loop over quadrature points
  for (int iQuadPt=0; iQuadPt < _numQuadPts; ++iQuadPt) {
    
    // Compute coordinates of quadrature point in cell
    // x = sum[i=0,n-1] (Ni * xi)
    // y = sum[i=0,n-1] (Ni * yi)
    for (int iBasis=0; iBasis < _numBasis; ++iBasis) {
      const double basis = _basisQuad[iQuadPt*_numBasis+iBasis];
      for (int iDim=0; iDim < _spaceDim; ++iDim)
	_quadPts[iQuadPt*_spaceDim+iDim] +=
	  basis * vertCoords[iBasis*_spaceDim+iDim];
    } // for
    
    // Compute Jacobian at quadrature point
    // J = [dx/dp dy/dp]
    // dx/dp = sum[i=0,n-1] (dNi/dp * xi)
    // dy/dp = sum[i=0,n-1] (dNi/dp * yi)
    for (int iBasis=0; iBasis < _numBasis; ++iBasis) {
      const double deriv = _basisDerivQuad[iQuadPt*_numBasis+iBasis];
      for (int iDim=0; iDim < _spaceDim; ++iDim)
	_jacobianQuad[iQuadPt*_spaceDim+iDim] += 
	  deriv * vertCoords[iBasis*_spaceDim+iDim];
    } // for

    // Compute determinant of Jacobian at quadrature point
    // |J| = sqrt(J transpose(J))
    double det = 0.0;
    for (int iDim=0; iDim < _spaceDim; ++iDim)
      det += _jacobianQuad[iQuadPt*_spaceDim+iDim] * 
	_jacobianQuad[iQuadPt*_spaceDim+iDim];
    det = sqrt(det);
    _checkJacobianDet(det);
    _jacobianDetQuad[iQuadPt] = det;

    // Compute inverse of Jacobian at quadrature point
    // Jinv = 1.0/[J]
    for (int iDim=0; iDim < _spaceDim; ++iDim)
      _jacobianInvQuad[iQuadPt*_spaceDim+iDim] = 
	1.0 / _jacobianQuad[iQuadPt*_spaceDim+iDim];
  } // for
} // computeGeometryQuad


// End of file 
