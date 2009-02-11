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
 * @file pylith/feassemble/Quadrature1D.hh
 *
 * @brief Quadrature for 1-D finite-elements.
 */

#if !defined(pylith_feassemble_quadrature1d_hh)
#define pylith_feassemble_quadrature1d_hh

#include "Quadrature.hh" // ISA Quadrature

template<typename mesh_type>
class pylith::feassemble::Quadrature1D : public Quadrature<mesh_type>
{ // Quadrature1D
  friend class TestQuadrature1D; // unit testing

// PUBLIC MEMBERS ///////////////////////////////////////////////////////
public :

  /// Constructor
  Quadrature1D(void);

  /// Destructor
  ~Quadrature1D(void);

  /// Create a copy of this object.
  Quadrature<mesh_type>* clone(void) const;

  /** Compute geometric quantities for a cell at quadrature points.
   *
   * @param mesh Finite-element mesh
   * @param coordinates Section containing vertex coordinates
   * @param cell Finite-element cell
   */
  void computeGeometry(const double* vertCoords,
                       const int coordDim,
                       const int cell);

// PROTECTED METHODS ////////////////////////////////////////////////////
protected :

  /** Copy constructor.
   *
   * @param q Quadrature to copy
   */
  Quadrature1D(const Quadrature1D& q);

// PRIVATE METHODS //////////////////////////////////////////////////////
private :

  const Quadrature1D& operator=(const Quadrature1D&); ///< Not implemented

}; // Quadrature1D

#include "Quadrature1D.icc" // inline methods
#include "Quadrature1D.cc" // template methods

#endif // pylith_feassemble_quadrature1d_hh

// End of file 
