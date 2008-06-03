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

#include "SlipTimeFn.hh"

/// Namespace for pylith package
namespace pylith {
  namespace faults {
    class BruneSlipFn;
    class TestBruneSlipFn; // unit testing
  } // faults
} // pylith

/// Namespace for spatialdata package
namespace spatialdata {
  namespace spatialdb {
    class SpatialDB;
  } // spatialdb
} // spatialdata

/// C++ implementation of Brune slip time function.
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

  /** Set spatial database for peak slip rate.
   *
   * @param db Spatial database
   */
  void dbPeakRate(spatialdata::spatialdb::SpatialDB* const db);

  /** Initialize slip time function.
   *
   * @param faultMesh Finite-element mesh of fault.
   * @param cs Coordinate system for mesh.
   * @param originTime Origin time for earthquake source.
   */
  void initialize(const ALE::Obj<Mesh>& faultMesh,
		  const spatialdata::geocoords::CoordSys* cs,
		  const double originTime =0.0);

  /** Get slip on fault surface at time t.
   *
   * @param slipField Slip field over fault surface.
   * @param t Time t.
   * @param faultMesh Mesh over fault surface.
   *
   * @returns Slip vector as left-lateral/reverse/normal.
   */
  void slip(const ALE::Obj<real_section_type>& slipField,
	    const double t,
	    const ALE::Obj<Mesh>& faultMesh);
  
  /** Get slip increment on fault surface between time t0 and t1.
   *
   * @param slipField Slip field over fault surface.
   * @param t0 Time t.
   * @param t1 Time t+dt.
   * @param faultMesh Mesh over fault surface.
   * 
   * @returns Increment in slip vector as left-lateral/reverse/normal.
   */
  void slipIncr(const ALE::Obj<real_section_type>& slipField,
		const double t0,
		const double t1,
		const ALE::Obj<Mesh>& faultMesh);


  /** Get final slip.
   *
   * @returns Final slip.
   */
  ALE::Obj<real_section_type> finalSlip(void);

  /** Get time when slip begins at each point.
   *
   * @returns Time when slip begins.
   */
  ALE::Obj<real_section_type> slipTime(void);

// NOT IMPLEMENTED //////////////////////////////////////////////////////
private :

  /// Not implemented
  BruneSlipFn(const BruneSlipFn& m);

  /// Not implemented
  const BruneSlipFn& operator=(const BruneSlipFn& f);

// PRIVATE METHODS //////////////////////////////////////////////////////
private :

  /** Compute slip using slip time function.
   *
   * @param t Time relative to slip starting time at point
   * @param finalSlip Final slip at point
   * @param peakRate Peak slip rate at point
   *
   * @returns Slip at point at time t
   */
  static
  double _slipFn(const double t,
		 const double finalSlip,
		 const double peakRate);

// PRIVATE MEMBERS //////////////////////////////////////////////////////
private :

  /// Parameters for Brune slip time function.
  /// Final slip (vector), peak slip rate (scalar), slip time (scalar).
  ALE::Obj<real_section_type> _parameters;

  /// Spatial database for final slip
  spatialdata::spatialdb::SpatialDB* _dbFinalSlip;

  /// Spatial database for slip time
  spatialdata::spatialdb::SpatialDB* _dbSlipTime;

   /// Spatial database for peak slip rate
  spatialdata::spatialdb::SpatialDB* _dbPeakRate;

  int _spaceDim; ///< Spatial dimension for slip field.

}; // class BruneSlipFn

#include "BruneSlipFn.icc" // inline methods

#endif // pylith_faults_bruneslipfn_hh


// End of file 
