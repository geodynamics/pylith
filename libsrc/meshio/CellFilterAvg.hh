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

// Include directives ---------------------------------------------------
#include "CellFilter.hh" // ISA CellFilter

// CellFilter -----------------------------------------------------------
template<typename mesh_type, typename field_type>
class pylith::meshio::CellFilterAvg : public CellFilter<mesh_type, field_type>
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
  CellFilter<mesh_type,field_type>* clone(void) const;

  /** Filter field over cells.
   *
   * @param fieldIn Field to filter.
   * @param label Label identifying cells.
   * @param labelId Value of label of cells to filter.
   *
   * @returns Averaged field.
   */
  const field_type&
  filter(const field_type& fieldIn,
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
  const CellFilterAvg& operator=(const CellFilterAvg&);

// PRIVATE MEMBERS //////////////////////////////////////////////////////
private :

  field_type* _fieldAvg; ///< Averaged cell field

}; // CellFilterAvg

#include "CellFilterAvg.cc" // template definitions

#endif // pylith_meshio_cellfilteravg_hh


// End of file 
