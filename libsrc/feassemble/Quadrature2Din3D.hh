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

/**
 * @file pylith/feassemble/Quadrature2Din3D.hh
 *
 * @brief Quadrature for 2-D finite-elements in 3-D space.
 */

#if !defined(pylith_feassemble_quadrature2din3d_hh)
#define pylith_feassemble_quadrature2din3d_hh

#include "Quadrature.hh"

namespace pylith {
  namespace feassemble {
    class Quadrature2Din3D;
    class TestQuadrature2Din3D;
  } // feassemble
} // pylith

class pylith::feassemble::Quadrature2Din3D : public Quadrature
{ // Quadrature2Din3D
  friend class TestQuadrature2Din3D; // unit testing

// PUBLIC MEMBERS ///////////////////////////////////////////////////////
public :

  /// Constructor
  Quadrature2Din3D(void);

  /// Destructor
  virtual ~Quadrature2Din3D(void);

  /// Create a copy of this object.
  virtual
  Quadrature* clone(void) const;

// PROTECTED METHODS ////////////////////////////////////////////////////
protected :

  /** Copy constructor.
   *
   * @param q Quadrature to copy
   */
  Quadrature2Din3D(const Quadrature2Din3D& q);

  /** Compute geometric quantities for a cell.
   *
   * @param coordinates Section containing vertex coordinates
   * @param cell Finite-element cell
   */
  void _computeGeometry(const ALE::Obj<ALE::Mesh::real_section_type>& coordinates,
			const ALE::Mesh::point_type& cell);

// PRIVATE METHODS //////////////////////////////////////////////////////
private :

  /// Not implemented
  const Quadrature2Din3D& operator=(const Quadrature2Din3D&);

}; // Quadrature2Din3D

#include "Quadrature2Din3D.icc" // inline methods

#endif // pylith_feassemble_quadrature2din3d_hh

// End of file 
