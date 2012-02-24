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

/** @file libsrc/faults/StepSlipFn.hh
 *
 * @brief C++ implementation of a step-function slip time function.
 */

#if !defined(pylith_faults_stepslipfn_hh)
#define pylith_faults_stepslipfn_hh

// Include directives ---------------------------------------------------
#include "SlipTimeFn.hh"

#include "pylith/utils/array.hh" // HASA double_array

// StepSlipFn -----------------------------------------------------------
/**
 * @brief C++ implementation of a step-function slip time function.
 *
 * Slip time function is a step function with slip beginning at time t0.
 *
 * Normalized slip = 1 if t >= t0, 0 otherwise
*/
class pylith::faults::StepSlipFn : public SlipTimeFn
{ // class StepSlipFn
  friend class TestStepSlipFn; // unit testing

// PUBLIC METHODS ///////////////////////////////////////////////////////
public :

  /// Default constructor.
  StepSlipFn(void);

  /// Destructor.
  ~StepSlipFn(void);

  /// Deallocate PETSc and local data structures.
  virtual
  void deallocate(void);
  
  /** Set spatial database for final slip.
   *
   * @param db Spatial database
   */
  void dbFinalSlip(spatialdata::spatialdb::SpatialDB* const db);

  /** Set spatial database for slip initiation time.
   *
   * @param db Spatial database
   */
  void dbSlipTime(spatialdata::spatialdb::SpatialDB* const db);

  /** Initialize slip time function.
   *
   * @param faultMesh Finite-element mesh of fault.
   * @param cs Coordinate system for mesh
   * @param normalizer Nondimensionalization of scales.
   * @param originTime Origin time for earthquake source.
   */
  void initialize(const topology::SubMesh& faultMesh,
		  const spatialdata::units::Nondimensional& normalizer,
		  const double originTime =0.0);

  /** Get slip on fault surface at time t.
   *
   * @param slipField Slip field over fault surface.
   * @param t Time t.
   *
   * @returns Slip vector as left-lateral/reverse/normal.
   */
  void slip(topology::Field<topology::SubMesh>* const slipField,
	    const double t);
  
  /** Get slip increment on fault surface between time t0 and t1.
   *
   * @param slipField Slip field over fault surface.
   * @param t0 Time t.
   * @param t1 Time t+dt.
   * 
   * @returns Increment in slip vector as left-lateral/reverse/normal.
   */
  void slipIncr(topology::Field<topology::SubMesh>* slipField,
		const double t0,
		const double t1);

  /** Get final slip.
   *
   * @returns Final slip.
   */
  const topology::Field<topology::SubMesh>& finalSlip(void);

  /** Get time when slip begins at each point.
   *
   * @returns Time when slip begins.
   */
  const topology::Field<topology::SubMesh>& slipTime(void);

// NOT IMPLEMENTED //////////////////////////////////////////////////////
private :

  StepSlipFn(const StepSlipFn&); ///< Not implemented.
  const StepSlipFn& operator=(const StepSlipFn&); ///< Not implemented

// PRIVATE MEMBERS //////////////////////////////////////////////////////
private :

  double _slipTimeVertex; ///< Slip time at a vertex.
  double_array _slipVertex; ///< Final slip at a vertex.

  /// Spatial database for final slip
  spatialdata::spatialdb::SpatialDB* _dbFinalSlip;

  /// Spatial database for slip time
  spatialdata::spatialdb::SpatialDB* _dbSlipTime;

}; // class StepSlipFn

#include "StepSlipFn.icc" // inline methods

#endif // pylith_faults_stepslipfn_hh


// End of file 
