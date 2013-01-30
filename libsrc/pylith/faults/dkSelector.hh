// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2012 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

/** @file libsrc/faults/dkSelector.hh
 *
 * @brief C++ implementation of the dynamic kinematic selector
 */

#if !defined(pylith_faults_dkselector_hh)
#define pylith_faults_dkselector_hh

// Include directives ---------------------------------------------------
#include "pylith/topology/topologyfwd.hh" // USES Fields<Field<SubMesh> >

#include "spatialdata/spatialdb/spatialdbfwd.hh" // USES SpatialDB

#include "pylith/utils/array.hh" // HASA scalar_array

// dkSelector -----------------------------------------------------------
/** @brief Dynamic-Kinematic Selector
 *
 * If value is under .5 the fault vertex has its slip value controlled by EQsrc
 *
 */
class pylith::faults::dkSelector 
{ // class dkSelector

// PUBLIC METHODS ///////////////////////////////////////////////////////
public :

  /// Default constructor.
  dkSelector(void);

  /// Destructor.
  ~dkSelector(void);

  /// Deallocate PETSc and local data structures.
  virtual
  void deallocate(void);
  
  /** Set spatial database for dkselector
   *
   * @param db dksel
   */
  void dbdksel(spatialdata::spatialdb::SpatialDB* const db);

  /** Initialize dkselector
   *
   * @param faultMesh Finite-element mesh of fault.
   * @param normalizer Nondimensionalization of scales.
   */
  void initialize(const topology::SubMesh& faultMesh,
		  const spatialdata::units::Nondimensional& normalizer);

  /** Get dk on fault surface (time will be implemented through this guy)
   *
   * @param dk DK selector field over fault surface
   *
   * @returns dk for the time
   */
  void dk(topology::Field<topology::SubMesh>* const dkField);
  
// NOT IMPLEMENTED //////////////////////////////////////////////////////
private :

  dkSelector(const dkSelector&); ///< Not implemented
  const dkSelector& operator=(const dkSelector&); ///< Not implemented

// PRIVATE MEMBERS //////////////////////////////////////////////////////
private :

  PylithScalar _dkselVertex; ///< Slip time at a vertex.

  /// Spatial database for slip rate.
  spatialdata::spatialdb::SpatialDB* _dbdksel;

}; // class dkSelector

#endif // pylith_faults_dkselector_hh


// End of file 
