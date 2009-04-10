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
 * @file pylith/meshio/VertexFilterVecNorm.hh
 *
 * @brief C++ object for computing vector norms for fields over
 * vertices when outputing finite-element data.
 */

#if !defined(pylith_meshio_cellfiltervecnorm_hh)
#define pylith_meshio_cellfiltervecnorm_hh

// Include directives ---------------------------------------------------
#include "VertexFilter.hh" // ISA VertexFilter

// VertexFilterVecNorm --------------------------------------------------
template<typename mesh_type>
class pylith::meshio::VertexFilterVecNorm : public VertexFilter<mesh_type>
{ // VertexFilterVecNorm

// PUBLIC METHODS ///////////////////////////////////////////////////////
public :

  /// Constructor
  VertexFilterVecNorm(void);

  /// Destructor
  ~VertexFilterVecNorm(void);

  /** Create copy of filter.
   *
   * @returns Copy of filter.
   */
  VertexFilter<mesh_type>* clone(void) const;

  /** Filter vertex field.
   *
   * @param fieldIn Field to filter.
   */
  const topology::Field<mesh_type>&
  filter(const topology::Field<mesh_type>& fieldIn);

// PROTECTED METHODS ////////////////////////////////////////////////////
protected :

  /** Copy constructor.
   *
   * @param f Filter to copy.
   * @returns Pointer to this.
   */
  VertexFilterVecNorm(const VertexFilterVecNorm& f);

// NOT IMPLEMENTED //////////////////////////////////////////////////////
private :

  /// Not implemented.
  const VertexFilterVecNorm& operator=(const VertexFilterVecNorm&);

// PRIVATE MEMBERS //////////////////////////////////////////////////////
private :

  topology::Field<mesh_type>* _fieldVecNorm; ///< Filtered vertex field

}; // VertexFilterVecNorm

#include "VertexFilterVecNorm.cc" // template definitions

#endif // pylith_meshio_cellfiltervecnorm_hh


// End of file 
