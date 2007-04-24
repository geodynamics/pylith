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

  /** Create copy of fault.
   *
   * @returns Copy of fault.
   */
  SlipTimeFn* clone(void) const;  

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
   * @param dbSlip Spatial database for slip.
   * @param dbSlipTime Spatial database for slip initiation time.
   * @param dbPeakRate Spatial database for peak slip rate.
   */
  void initialize(const ALE::Obj<Mesh>& mesh,
		  const spatialdata::geocoords::CoordSys* cs,
		  const std::set<Mesh::point_type>& vertices);

  /** Get slip on fault surface at time t.
   *
   * @param t Time t.
   * @param vertices Vertices on fault surface.
   */
  const ALE::Obj<real_section_type>& slip(const double t,
				 const std::set<Mesh::point_type>& vertices);

// PROTECTED METHODS ////////////////////////////////////////////////////
protected :

  /** Copy constructor.
   *
   * @param f BruneSlipFn to copy
   */
  BruneSlipFn(const BruneSlipFn& m);

// NOT IMPLEMENTED //////////////////////////////////////////////////////
private :

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
  double _slip(const double t,
	       const double finalSlip,
	       const double peakRate);

// PRIVATE MEMBERS //////////////////////////////////////////////////////
private :

  ALE::Obj<real_section_type> _slipField; ///< Slip field on fault surface

  /// Spatial database for final slip
  spatialdata::spatialdb::SpatialDB* _dbFinalSlip;

  /// Spatial database for slip time
  spatialdata::spatialdb::SpatialDB* _dbSlipTime;

   /// Spatial database for peak slip rate
  spatialdata::spatialdb::SpatialDB* _dbPeakRate;

}; // class BruneSlipFn

#include "BruneSlipFn.icc" // inline methods

#endif // pylith_faults_bruneslipfn_hh


// End of file 
