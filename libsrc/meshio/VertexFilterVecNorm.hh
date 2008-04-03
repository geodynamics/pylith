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

#include "VertexFilter.hh" // ISA VertexFilter

namespace pylith {
  namespace meshio {
    class VertexFilterVecNorm;
  } // meshio
} // pylith

class pylith::meshio::VertexFilterVecNorm : public VertexFilter
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
  VertexFilter* clone(void) const;

  /** Filter field. Field type of filtered field is returned via an argument.
   *
   * @param fieldType Field type of filtered field.
   * @param fieldIn Field to filter.
   * @param mesh PETSc mesh.
   */
  const ALE::Obj<real_section_type>&
  filter(VectorFieldEnum* fieldType,
	 const ALE::Obj<real_section_type>& fieldIn,
	 const ALE::Obj<Mesh>& mesh);

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
  const VertexFilter& operator=(const VertexFilter&);

// PRIVATE MEMBERS //////////////////////////////////////////////////////
private :

  ALE::Obj<real_section_type> _fieldVecNorm; ///< Filtered vertex field

}; // VertexFilterVecNorm

#endif // pylith_meshio_cellfiltervecnorm_hh


// End of file 
