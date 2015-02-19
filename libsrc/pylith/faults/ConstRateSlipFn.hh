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
// Copyright (c) 2010-2015 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

/** @file libsrc/faults/ConstRateSlipFn.hh
 *
 * @brief C++ implementation of a constant slip rate slip time function.
 */

#if !defined(pylith_faults_constrateslipfn_hh)
#define pylith_faults_constrateslipfn_hh

// Include directives ---------------------------------------------------
#include "SlipTimeFn.hh"

#include "pylith/topology/topologyfwd.hh" // USES Fields<Field<Mesh> >

#include "spatialdata/spatialdb/spatialdbfwd.hh" // USES SpatialDB

#include "pylith/utils/array.hh" // HASA scalar_array

// ConstRateTimeFn ------------------------------------------------------
/** @brief Constant slip rate slip-time function.
 *
 * Slip time function follows the integral of constant slip rate slip
 * time function.
 *
 * Normalized slip = sliprate * (t - t0)
 */
class pylith::faults::ConstRateSlipFn : public SlipTimeFn
{ // class ConstRateSlipFn
  friend class TestConstRateSlipFn; // unit testing

// PUBLIC METHODS ///////////////////////////////////////////////////////
public :

  /// Default constructor.
  ConstRateSlipFn(void);

  /// Destructor.
  ~ConstRateSlipFn(void);

  /// Deallocate PETSc and local data structures.
  virtual
  void deallocate(void);
  
  /** Set spatial database for slip rate.
   *
   * @param db Spatial database
   */
  void dbSlipRate(spatialdata::spatialdb::SpatialDB* const db);

  /** Set spatial database for slip initiation time.
   *
   * @param db Spatial database
   */
  void dbSlipTime(spatialdata::spatialdb::SpatialDB* const db);

  /** Initialize slip time function.
   *
   * @param faultMesh Finite-element mesh of fault.
   * @param normalizer Nondimensionalization of scales.
   * @param originTime Origin time for earthquake source.
   */
  void initialize(const topology::Mesh& faultMesh,
		  const spatialdata::units::Nondimensional& normalizer,
		  const PylithScalar originTime =0.0);

  /** Get slip on fault surface at time t.
   *
   * @param slipField Slip field over fault surface.
   * @param t Time t.
   *
   * @returns Slip vector as left-lateral/reverse/normal.
   */
  void slip(topology::Field* const slipField,
	    const PylithScalar t);
  
  /** Get final slip.
   *
   * @returns Final slip.
   */
  const topology::Field& finalSlip(void);

  /** Get time when slip begins at each point.
   *
   * @returns Time when slip begins.
   */
  const topology::Field& slipTime(void);

// NOT IMPLEMENTED //////////////////////////////////////////////////////
private :

  ConstRateSlipFn(const ConstRateSlipFn&); ///< Not implemented
  const ConstRateSlipFn& operator=(const ConstRateSlipFn&); ///< Not implemented

// PRIVATE MEMBERS //////////////////////////////////////////////////////
private :

  PylithScalar _slipTimeVertex; ///< Slip time at a vertex.
  scalar_array _slipRateVertex; ///< Slip rate at a vertex.

  /// Spatial database for slip rate.
  spatialdata::spatialdb::SpatialDB* _dbSlipRate;

  /// Spatial database for slip time.
  spatialdata::spatialdb::SpatialDB* _dbSlipTime;

}; // class ConstRateSlipFn

#include "ConstRateSlipFn.icc" // inline methods

#endif // pylith_faults_constrateslipfn_hh


// End of file 
