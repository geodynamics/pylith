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
#include "pylith/utils/vectorfields.hh" // USES VectorFieldEnum

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

  /** Filter field. Field type of filtered field is returned via an argument.
   *
   * @param fieldType Field type of filtered field.
   * @param fieldIn Field to filter.
   * @param mesh PETSc mesh.
   * @param label Label identifying cells.
   * @param Value of label of cells to filter.
   *
   * @returns Averaged field.
   */
  virtual
  const ALE::Obj<real_section_type>&
  filter(VectorFieldEnum* fieldType,
	 const ALE::Obj<real_section_type>& fieldIn,
	 const ALE::Obj<ALE::Mesh>& mesh,
	 const char* label =0,
	 const int labelId =0) = 0;

// PROTECTED METHODS ////////////////////////////////////////////////////
protected :

  /** Copy constructor.
   *
   * @param f Filter to copy.
   * @returns Pointer to this.
   */
  CellFilter(const CellFilter& f);

// NOT IMPLEMENTED //////////////////////////////////////////////////////
private :

  /// Not implemented.
  const CellFilter& operator=(const CellFilter&);

// PROTECTED MEMBERS ////////////////////////////////////////////////////
protected :

  feassemble::Quadrature* _quadrature; ///< Quadrature associated with cells.

}; // CellFilter

#endif // pylith_meshio_cellfilter_hh


// End of file 
