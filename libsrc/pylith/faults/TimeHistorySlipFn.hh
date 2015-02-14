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

/** @file libsrc/faults/TimeHistorySlipFn.hh
 *
 * @brief C++ implementation of a slip time function with a
 * user-specified time history.
 */

#if !defined(pylith_faults_timehistoryslipfn_hh)
#define pylith_faults_timehistoryslipfn_hh

// Include directives ---------------------------------------------------
#include "SlipTimeFn.hh"

#include "spatialdata/spatialdb/spatialdbfwd.hh" // USES SpatialDB

#include "pylith/utils/array.hh" // HASA scalar_array

// TimeHistorySlipFn -----------------------------------------------------------
/**
 * @brief C++ implementation of a slip time function with a
 * user-specified time history.
 *
 * User-specified slip time function with spatial variable amplitude
 * and starting time t0.
*/
class pylith::faults::TimeHistorySlipFn : public SlipTimeFn
{ // class TimeHistorySlipFn
  friend class TestTimeHistorySlipFn; // unit testing

// PUBLIC METHODS ///////////////////////////////////////////////////////
public :

  /// Default constructor.
  TimeHistorySlipFn(void);

  /// Destructor.
  ~TimeHistorySlipFn(void);

  /// Deallocate PETSc and local data structures.
  virtual
  void deallocate(void);
  
  /** Set spatial database for amplitude of slip time history.
   *
   * @param db Spatial database
   */
  void dbAmplitude(spatialdata::spatialdb::SpatialDB* const db);

  /** Set spatial database for slip initiation time.
   *
   * @param db Spatial database
   */
  void dbSlipTime(spatialdata::spatialdb::SpatialDB* const db);

  /** Set time history.
   *
   * @param th Time history.
   */
  void dbTimeHistory(spatialdata::spatialdb::TimeHistory* const th);

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

  TimeHistorySlipFn(const TimeHistorySlipFn&); ///< Not implemented.
  const TimeHistorySlipFn& operator=(const TimeHistorySlipFn&); ///< Not implemented

// PRIVATE MEMBERS //////////////////////////////////////////////////////
private :

  PylithScalar _slipTimeVertex; ///< Slip time at a vertex.
  PylithScalar _timeScale; ///< Time scale.
  scalar_array _slipVertex; ///< Final slip at a vertex.

  /// Spatial database for amplitude of slip time history.
  spatialdata::spatialdb::SpatialDB* _dbAmplitude;

  /// Spatial database for slip time.
  spatialdata::spatialdb::SpatialDB* _dbSlipTime;

  /// Time history database.
  spatialdata::spatialdb::TimeHistory* _dbTimeHistory;

}; // class TimeHistorySlipFn

#include "TimeHistorySlipFn.icc" // inline methods

#endif // pylith_faults_timehistoryslipfn_hh


// End of file 
