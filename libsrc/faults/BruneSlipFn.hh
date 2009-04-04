// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
//                           Brad T. Aagaard
//                        U.S. Geological Survey
//
// {LicenseText}
//
// ----------------------------------------------------------------------
//

/** @file libsrc/faults/BruneSlipFn.hh
 *
 * @brief C++ implementation of Brune slip time function.
 *
 * Slip time function follows the integral of Brune's (1907) far-field
 * time function.
 *
 * Normalize slip = 1 - exp(-t/tau)(1 + t/tau),
 * where tau = finalSlip / (exp(1.0) * peakRate)
 */

#if !defined(pylith_faults_bruneslipfn_hh)
#define pylith_faults_bruneslipfn_hh

// Include directives ---------------------------------------------------
#include "SlipTimeFn.hh"

#include "pylith/topology/topologyfwd.hh" // USES Fields<Field<SubMesh> >

#include "pylith/utils/array.hh" // HASA double_array

// SlipTimeFn -----------------------------------------------------------
class pylith::faults::BruneSlipFn : public SlipTimeFn
{ // class BruneSlipFn
  friend class TestBruneSlipFn; // unit testing

// PUBLIC METHODS ///////////////////////////////////////////////////////
public :

  /// Default constructor.
  BruneSlipFn(void);

  /// Destructor.
  virtual
  ~BruneSlipFn(void);

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

  /** Set spatial database for rise time (0 -> 0.95 final slip).
   *
   * @param db Spatial database
   */
  void dbRiseTime(spatialdata::spatialdb::SpatialDB* const db);

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

  BruneSlipFn(const BruneSlipFn&); ///< Not implemented
  const BruneSlipFn& operator=(const BruneSlipFn&); ///< Not implemented

// PRIVATE METHODS //////////////////////////////////////////////////////
private :

  /** Compute slip using slip time function.
   *
   * @param t Time relative to slip starting time at point.
   * @param finalSlip Final slip at point.
   * @param riseTime Peak slip rate at point.
   *
   * @returns Slip at point at time t
   */
  static
  double _slipFn(const double t,
		 const double finalSlip,
		 const double riseTime);

// PRIVATE MEMBERS //////////////////////////////////////////////////////
private :

  double _slipTimeVertex; ///< Slip time at a vertex.
  double _riseTimeVertex; ///< Rise time at a vertex.
  double_array _finalSlipVertex; ///< Final slip at a vertex.

  /// Parameters for Brune slip time function, final slip (vector),
  /// rise time (scalar), slip time (scalar).
  topology::Fields<topology::Field<topology::SubMesh> >* _parameters;

  /// Spatial database for final slip.
  spatialdata::spatialdb::SpatialDB* _dbFinalSlip;

  /// Spatial database for slip time.
  spatialdata::spatialdb::SpatialDB* _dbSlipTime;

   /// Spatial database for rise time (0 -> 0.95 final slip).
  spatialdata::spatialdb::SpatialDB* _dbRiseTime;

}; // class BruneSlipFn

#include "BruneSlipFn.icc" // inline methods

#endif // pylith_faults_bruneslipfn_hh


// End of file 
