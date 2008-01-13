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
 * @file pylith/meshio/OutputFilter.hh
 *
 * @brief C++ object for filtering finite-element output fields.
 */

#if !defined(pylith_meshio_outputfilter_hh)
#define pylith_meshio_outputfilter_hh

#include "pylith/utils/sievetypes.hh" // USES ALE::Mesh, real_section_type

namespace pylith {
  namespace meshio {
    class OutputFilter;
  } // meshio

} // pylith

class pylith::meshio::OutputFilter
{ // OutputFilter

// PUBLIC ENUMS /////////////////////////////////////////////////////////
public :

  /// Type of filter.
  enum FilterEnum { 
    VERTEX_FILTER, ///< Filter vertex field.
    CELL_FILTER ///< Filter cell field.
  };

// PUBLIC METHODS ///////////////////////////////////////////////////////
public :

  /// Constructor
  OutputFilter(void);

  /// Destructor
  ~OutputFilter(void);

  /** Create copy of filter.
   *
   * @returns Copy of filter.
   */
  virtual
  OutputFilter* clone(void) const = 0;

  /** Get filter type.
   *
   * @returns Type of filter.
   */
  FilterEnum filterType(void) const;

  /** Filter field.
   *
   * @param fieldIn Field to filter.
   * @param mesh PETSc mesh.
   */
  virtual
  const ALE::Obj<real_section_type>&
  filter(const ALE::Obj<real_section_type>& fieldIn,
	 const ALE::Obj<ALE::Mesh>& mesh) = 0;

// PROTECTED METHODS ////////////////////////////////////////////////////
protected :

  /** Copy constructor.
   *
   * @param f Filter to copy.
   * @returns Pointer to this.
   */
  OutputFilter(const OutputFilter& f);

  /** operator=.
  *
  * @param f Filter to copy.
  * @returns Copy of filter.
  */
  const OutputFilter& operator=(const OutputFilter& f);

// PRIVATE MEMBERS //////////////////////////////////////////////////////
private :

  FilterEnum _filterType; ///< Type of filter.

}; // OutputFilter

#endif // pylith_meshio_outputfilter_hh


// End of file 
