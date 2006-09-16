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

#if !defined(pylith_feassemble_quadrature1d_hh)
#define pylith_feassemble_quadrature1d_hh

#include "Quadrature.hh"

namespace pylith {
  namespace feassemble {
    class Quadrature1D;
  } // feassemble
} // pylith

class pylith::feassemble::Quadrature1D
{ // Quadrature1D
  
// PUBLIC MEMBERS ///////////////////////////////////////////////////////
public :

  /// Constructor
  Quadrature1D(void);

  /// Destructor
  ~Quadrature1D(void);

  /** Compute geometry of element object.
   *
   */
  void compute(const ALE::Obj<ALE::Mesh::section_type>& coordinates,
	       const ALE::Mesh::point_type& cell,
	       value_type* pV,
	       value_type* pJacobian,
	       valye_type* pJacobianInv,
	       value_type* pJacobianDet);

}; // Quadrature1D

#endif // pylith_feassemble_quadrature1d_hh

// End of file 
