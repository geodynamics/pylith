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

#include "ElemGeometry1D.hh" // implementation of class methods

// ----------------------------------------------------------------------
// Constructor
pylith::feassemble::ElemGeometry1D::ElemGeometry1D(void)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor
pylith::feassemble::ElemGeometry1D::~ElemGeometry1D(void)
{ // destructor
} // destructor
  
// ----------------------------------------------------------------------
/** Compute geometry of element object.
 *
 */
void compute(const Obj<section_type>& coordinates,
	     const point_type& cell,
	     value_type* pV,
	     value_type* pJacobian,
	     valye_type* pJacobianInv,
	     value_type& jacobianDet)
{ // compute
  assert(0 != pV);
  assert(0 != pJacobian);
  assert(0 != pJacobianInv);

  const int ndims = 1;
  const topology_type::patch_type patch  = 0;
  const value_type* coords = coordinates->restrict(patch, cell);

  // Matt: What is pV???
  // Are they the coordinates of the vertices of the cell?
  pV[0] = coords[0];

  


  value_type* invDet;





  
  for(int d = 0; d < _dim; d++) {
    for(int f = 0; f < _dim; f++) {
      J[d*_dim+f] = 0.5*(coords[(f+1)*_dim+d] - coords[0*_dim+d]);
    }
  }
  if (_dim == 1) {
    detJ = J[0];
  } else if (_dim == 2) {
    detJ = J[0]*J[3] - J[1]*J[2];
  } else if (_dim == 3) {
    detJ = J[0*3+0]*(J[1*3+1]*J[2*3+2] - J[1*3+2]*J[2*3+1]) +
      J[0*3+1]*(J[1*3+2]*J[2*3+0] - J[1*3+0]*J[2*3+2]) +
      J[0*3+2]*(J[1*3+0]*J[2*3+1] - J[1*3+1]*J[2*3+0]);
  }
  if (invJ) {
    invDet = 1.0/detJ;
    if (_dim == 2) {
      invJ[0] =  invDet*J[3];
      invJ[1] = -invDet*J[1];
      invJ[2] = -invDet*J[2];
      invJ[3] =  invDet*J[0];
    } else if (_dim == 3) {
      // FIX: This may be wrong
      invJ[0*3+0] = invDet*(J[1*3+1]*J[2*3+2] - J[1*3+2]*J[2*3+1]);
      invJ[0*3+1] = invDet*(J[1*3+2]*J[2*3+0] - J[1*3+0]*J[2*3+2]);
      invJ[0*3+2] = invDet*(J[1*3+0]*J[2*3+1] - J[1*3+1]*J[2*3+0]);
      invJ[1*3+0] = invDet*(J[0*3+1]*J[2*3+2] - J[0*3+2]*J[2*3+1]);
      invJ[1*3+1] = invDet*(J[0*3+2]*J[2*3+0] - J[0*3+0]*J[2*3+2]);
      invJ[1*3+2] = invDet*(J[0*3+0]*J[2*3+1] - J[0*3+1]*J[2*3+0]);
      invJ[2*3+0] = invDet*(J[0*3+1]*J[1*3+2] - J[0*3+2]*J[1*3+1]);
      invJ[2*3+1] = invDet*(J[0*3+2]*J[1*3+0] - J[0*3+0]*J[1*3+2]);
      invJ[2*3+2] = invDet*(J[0*3+0]*J[1*3+1] - J[0*3+1]*J[1*3+0]);
    }
  }
  
} // compute

// End of file 
