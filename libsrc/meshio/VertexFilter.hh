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

#include "pylith/utils/sievetypes.hh" // USES PETSc Mesh, real_section_type
#include "pylith/utils/vectorfields.hh" // USES VectorFieldEnum

namespace pylith {
  namespace meshio {
    class VertexFilter;
  } // meshio

} // pylith

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

  /** Filter field. Field type of filtered field is returned via an argument.
   *
   * @param fieldType Field type of filtered field.
   * @param fieldIn Field to filter.
   * @param mesh PETSc mesh.
   */
  virtual
  const ALE::Obj<real_section_type>&
  filter(VectorFieldEnum* fieldType,
	 const ALE::Obj<real_section_type>& fieldIn,
	 const ALE::Obj<Mesh>& mesh) = 0;

// PROTECTED METHODS ////////////////////////////////////////////////////
protected :

  /** Copy constructor.
   *
   * @param f Filter to copy.
   * @returns Pointer to this.
   */
  VertexFilter(const VertexFilter& f);

  /** operator=.
  *
  * @param f Filter to copy.
  * @returns Copy of filter.
  */
  const VertexFilter& operator=(const VertexFilter& f);

}; // VertexFilter

#endif // pylith_meshio_vertexfilter_hh


// End of file 
