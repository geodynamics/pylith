// -*- C++ -*-
//
// ======================================================================
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2015 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

/**
 * @file libsrc/meshio/CellFilterAvg.hh
 *
 * @brief C++ object for averaging cell fields over quadrature points
 * when outputing finite-element data.
 */

#if !defined(pylith_meshio_cellfilteravg_hh)
#define pylith_meshio_cellfilteravg_hh

// Include directives ---------------------------------------------------
#include "CellFilter.hh" // ISA CellFilter

// CellFilter -----------------------------------------------------------
/** @brief C++ object for averaging cell fields over quadrature points
 * when outputing finite-element data.
 */
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

  /// Deallocate PETSc and local data structures.
  void deallocate(void);
  
  /** Get averaged field buffer.
   *
   * @returns Field associated with averaged values.
   */
  const topology::Field* fieldAvg(void) const;
  
  /** Filter field over cells.
   *
   * @param fieldIn Field to filter.
   * @param label Label identifying cells.
   * @param labelId Value of label of cells to filter.
   *
   * @returns Averaged field.
   */
  topology::Field&
  filter(const topology::Field& fieldIn,
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

  topology::Field* _fieldAvg; ///< Averaged cell field

}; // CellFilterAvg


#endif // pylith_meshio_cellfilteravg_hh


// End of file 
