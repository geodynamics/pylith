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
 * @file libsrc/meshio/CellFilter.hh
 *
 * @brief C++ object for filtering cell fields when outputing
 * finite-element data.
 */

#if !defined(pylith_meshio_cellfilter_hh)
#define pylith_meshio_cellfilter_hh

// Include directives ---------------------------------------------------
#include "meshiofwd.hh" // forward declarations

#include "pylith/feassemble/feassemblefwd.hh" // HOLDSA Quadrature
#include "pylith/topology/topologyfwd.hh" // USES Field

// CellFilter -----------------------------------------------------------
/** @brief C++ object for filtering cell fields when outputing
 * finite-element data.
 */
class pylith::meshio::CellFilter
{ // CellFilter

// PUBLIC METHODS ///////////////////////////////////////////////////////
public :

  /// Constructor
  CellFilter(void);

  /// Destructor
  virtual
  ~CellFilter(void);

  /** Create copy of filter.
   *
   * @returns Copy of filter.
   */
  virtual
  CellFilter* clone(void) const = 0;

  /// Deallocate PETSc and local data structures.
  virtual
  void deallocate(void);
  
  /** Set quadrature associated with cells.
   *
   * @param q Quadrature for cells.
   */
  void quadrature(const feassemble::Quadrature* q);

  /** Filter field. Field type of filtered field is returned via an argument.
   *
   * @param fieldIn Field to filter.
   * @param label Value of label of cells to filter.
   * @param labelId Id associated with label of cells to filter.
   *
   * @returns Averaged field.
   */
  virtual
  topology::Field&
  filter(const topology::Field& fieldIn,
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

  /// Quadrature associated with cells.
  feassemble::Quadrature* _quadrature;

  topology::StratumIS* _cellsIS; ///< Index set of cells.

}; // CellFilter

#endif // pylith_meshio_cellfilter_hh


// End of file 
