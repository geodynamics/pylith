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
 * @file pylith/meshio/CellFilterAvg.hh
 *
 * @brief C++ object for averaging cell fields over quadrature points
 * when outputing finite-element data.
 */

#if !defined(pylith_meshio_cellfilteravg_hh)
#define pylith_meshio_cellfilteravg_hh

#include "CellFilter.hh" // ISA CellFilter

namespace pylith {
  namespace meshio {
    class CellFilterAvg;
  } // meshio
} // pylith

class pylith::meshio::CellFilterAvg : public CellFilter
{ // CellFilterAvg

// PUBLIC METHODS ///////////////////////////////////////////////////////
public :

  /// Constructor
  CellFilterAvg(void);

  /// Destructor
  ~CellFilterAvg(void);

  /** Create copy of filter.
   *
   * @returns Copy of filter.
   */
  CellFilter* clone(void) const;

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
  const ALE::Obj<real_section_type>&
  filter(VectorFieldEnum* fieldType,
	 const ALE::Obj<real_section_type>& fieldIn,
	 const ALE::Obj<ALE::Mesh>& mesh,
	 const char* label =0,
	 const int labelId =0);

// PROTECTED METHODS ////////////////////////////////////////////////////
protected :

  /** Copy constructor.
   *
   * @param f Filter to copy.
   * @returns Pointer to this.
   */
  CellFilterAvg(const CellFilterAvg& f);

// NOT IMPLEMENTED //////////////////////////////////////////////////////
private :

  /// Not implemented.
  const CellFilter& operator=(const CellFilter&);

// PRIVATE MEMBERS //////////////////////////////////////////////////////
private :

  ALE::Obj<real_section_type> _fieldAvg; ///< Averaged cell field

}; // CellFilterAvg

#endif // pylith_meshio_cellfilteravg_hh


// End of file 
