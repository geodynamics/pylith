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

#if !defined(pylith_feassemble_elemgeometry1d_hh)
#define pylith_feassemble_elemgeometry1d_hh

#include "ElemGeometry.hh"

namespace pylith {
  namespace feassemble {
    class ElemGeometry1D;
  } // feassemble
} // pylith

class pylith::feassemble::ElemGeometry1D
{ // ElemGeometry1D
  
// PUBLIC MEMBERS ///////////////////////////////////////////////////////
public :

  /// Constructor
  ElemGeometry1D(void);

  /// Destructor
  ~ElemGeometry1D(void);

  /** Compute geometry of element object.
   *
   */
  void compute(const Obj<section_type>& coordinates,
	       const point_type& cell,
	       value_type* pV,
	       value_type* pJacobian,
	       valye_type* pJacobianInv,
	       value_type* pJacobianDet);

}; // ElemGeometry1D

#endif // pylith_feassemble_elemgeometry1d_hh

// End of file 
