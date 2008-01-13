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
 * @file pylith/meshio/CellFilter.hh
 *
 * @brief C++ object for filtering cell fields when outputing
 * finite-element data.
 */

#if !defined(pylith_meshio_cellfilter_hh)
#define pylith_meshio_cellfilter_hh

#include "pylith/utils/sievetypes.hh" // USES ALE::Mesh, real_section_type

namespace pylith {
  namespace meshio {
    class CellFilter;
  } // meshio

  namespace feassemble {
    class Quadrature;
  } // meshio  
} // pylith

class pylith::meshio::CellFilter
{ // CellFilter

// PUBLIC METHODS ///////////////////////////////////////////////////////
public :

  /// Constructor
  CellFilter(void);

  /// Destructor
  ~CellFilter(void);

  /** Create copy of filter.
   *
   * @returns Copy of filter.
   */
  virtual
  CellFilter* clone(void) const = 0;

  /** Set quadrature associated with cells.
   *
   * @param q Quadrature for cells.
   */
  void quadrature(const feassemble::Quadrature* q);

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
  CellFilter(const CellFilter& f);

  /** operator=.
  *
  * @param f Filter to copy.
  * @returns Copy of filter.
  */
  const CellFilter& operator=(const CellFilter& f);

// PROTECTED MEMBERS ////////////////////////////////////////////////////
protected :

  feassemble::Quadrature* _quadrature; ///< Quadrature associated with cells.

}; // CellFilter

#endif // pylith_meshio_cellfilter_hh


// End of file 
