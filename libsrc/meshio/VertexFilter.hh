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
 * @file pylith/meshio/VertexFilter.hh
 *
 * @brief C++ object for filtering vertex fields when outputing
 * finite-element data.
 */

#if !defined(pylith_meshio_vertexfilter_hh)
#define pylith_meshio_vertexfilter_hh

// Include directives ---------------------------------------------------
#include "meshiofwd.hh" // forward declarations

#include "pylith/topology/topologyfwd.hh" // USES Field

// VertexFilter ---------------------------------------------------------
template<typename mesh_type>
class pylith::meshio::VertexFilter
{ // VertexFilter

// PUBLIC METHODS ///////////////////////////////////////////////////////
public :

  /// Constructor
  VertexFilter(void);

  /// Destructor
  ~VertexFilter(void);

  /** Create copy of filter.
   *
   * @returns Copy of filter.
   */
  virtual
  VertexFilter* clone(void) const = 0;

  /** Filter field over vertices of a mesh.
   *
   * @param fieldIn Field to filter.
   */
  virtual
  const topology::Field<mesh_type>&
  filter(const topology::Field<mesh_type>& fieldIn) = 0;

// PROTECTED METHODS ////////////////////////////////////////////////////
protected :

  /** Copy constructor.
   *
   * @param f Filter to copy.
   * @returns Pointer to this.
   */
  VertexFilter(const VertexFilter& f);

private :
  /** operator=.
  *
  * @param f Filter to copy.
  * @returns Copy of filter.
  */
  const VertexFilter& operator=(const VertexFilter& f);

}; // VertexFilter

#include "VertexFilter.cc" // template definitions

#endif // pylith_meshio_vertexfilter_hh


// End of file 
