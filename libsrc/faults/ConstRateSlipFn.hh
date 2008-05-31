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

/** @file libsrc/faults/ConstRateSlipFn.hh
 *
 * @brief C++ implementation of a constant slip rate slip time function.
 *
 * Slip time function follows the integral of constant slip rate slip
 * time function.
 *
 * Normalized slip = sliprate * (t - t0)
 */

#if !defined(pylith_faults_constrateslipfn_hh)
#define pylith_faults_constrateslipfn_hh

#include "SlipTimeFn.hh"

/// Namespace for pylith package
namespace pylith {
  namespace faults {
    class ConstRateSlipFn;
    class TestConstRateSlipFn; // unit testing
  } // faults
} // pylith

/// Namespace for spatialdata package
namespace spatialdata {
  namespace spatialdb {
    class SpatialDB;
  } // spatialdb
} // spatialdata

/// C++ implementation of ConstRate slip time function.
class pylith::faults::ConstRateSlipFn : public SlipTimeFn
{ // class ConstRateSlipFn
  friend class TestConstRateSlipFn; // unit testing

// PUBLIC METHODS ///////////////////////////////////////////////////////
public :

  /// Default constructor.
  ConstRateSlipFn(void);

  /// Destructor.
  virtual
  ~ConstRateSlipFn(void);

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
   * @param cs Coordinate system for mesh.
   * @param originTime Origin time for earthquake source.
   */
  void initialize(const ALE::Obj<Mesh>& faultMesh,
		  const spatialdata::geocoords::CoordSys* cs,
		  const double originTime =0.0);

  /** Get slip on fault surface at time t.
   *
   * @param t Time t.
   * @param faultMesh Mesh over fault surface.
   */
  const ALE::Obj<real_section_type>& slip(const double t,
					  const ALE::Obj<Mesh>& faultMesh);

  /** Get slip increment on fault surface between time t0 and t1.
   *
   * @param t0 Time t.
   * @param t1 Time t+dt.
   * @param faultMesh Mesh over fault surface.
   */
  const ALE::Obj<real_section_type>& slipIncr(const double t0,
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
  ConstRateSlipFn(const ConstRateSlipFn& m);

  /// Not implemented
  const ConstRateSlipFn& operator=(const ConstRateSlipFn& f);

// PRIVATE MEMBERS //////////////////////////////////////////////////////
private :

  /// Parameters for ConstRate slip time function.
  /// Final slip (vector), peak slip rate (scalar), slip time (scalar).
  ALE::Obj<real_section_type> _parameters;

  ALE::Obj<real_section_type> _slip; ///< Slip field on fault surface

  /// Spatial database for final slip
  spatialdata::spatialdb::SpatialDB* _dbSlipRate;

  /// Spatial database for slip time
  spatialdata::spatialdb::SpatialDB* _dbSlipTime;

  int _spaceDim; ///< Spatial dimension for slip field.

}; // class ConstRateSlipFn

#include "ConstRateSlipFn.icc" // inline methods

#endif // pylith_faults_constrateslipfn_hh


// End of file 
