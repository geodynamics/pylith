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

/** @file libsrc/faults/LiuCosSlipFn.hh
 *
 * @brief C++ implementation of LiuCos slip time function.
 *
 * Sine/cosine slip time function from Liu, Archuleta, and Hartzell,
 * BSSA, 2006 (doi:10.1785/0120060036) which has a rapid rise and then
 * a gradual falloff with a finite duration.
 */

#if !defined(pylith_faults_liucosslipfn_hh)
#define pylith_faults_liucosslipfn_hh

#include "SlipTimeFn.hh"

/// Namespace for pylith package
namespace pylith {
  namespace faults {
    class LiuCosSlipFn;
    class TestLiuCosSlipFn; // unit testing
  } // faults
} // pylith

/// Namespace for spatialdata package
namespace spatialdata {
  namespace spatialdb {
    class SpatialDB;
  } // spatialdb
  namespace units {
    class Nondimensional;
  } // units
} // spatialdata

/// C++ implementation of LiuCos slip time function.
class pylith::faults::LiuCosSlipFn : public SlipTimeFn
{ // class LiuCosSlipFn
  friend class TestLiuCosSlipFn; // unit testing

// PUBLIC METHODS ///////////////////////////////////////////////////////
public :

  /// Default constructor.
  LiuCosSlipFn(void);

  /// Destructor.
  virtual
  ~LiuCosSlipFn(void);

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

  /** Set spatial database for rise time. The rise time is the time it
   * takes for the slip to increase from 0.0 to 0.95 of the final
   * value.
   *
   * @param db Spatial database
   */
  void dbRiseTime(spatialdata::spatialdb::SpatialDB* const db);

  /** Initialize slip time function.
   *
   * @param faultMesh Finite-element mesh of fault.
   * @param cs Coordinate system for mesh.
   * @param originTime Origin time for earthquake source.
   * @param normalizer Nondimensionalization of scales.
   */
  void initialize(const ALE::Obj<Mesh>& faultMesh,
		  const spatialdata::geocoords::CoordSys* cs,
		  const spatialdata::units::Nondimensional& normalizer,
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
  LiuCosSlipFn(const LiuCosSlipFn& m);

  /// Not implemented
  const LiuCosSlipFn& operator=(const LiuCosSlipFn& f);

// PRIVATE METHODS //////////////////////////////////////////////////////
private :

  /** Compute slip using slip time function.
   *
   * @param t Time relative to slip starting time at point
   * @param finalSlip Final slip at point
   * @param riseTime Rise time (t95) at point
   *
   * @returns Slip at point at time t
   */
  static
  double _slipFn(const double t,
		 const double finalSlip,
		 const double riseTime);

// PRIVATE MEMBERS //////////////////////////////////////////////////////
private :

  /// Parameters for LiuCos slip time function.
  /// Final slip (vector), peak slip rate (scalar), slip time (scalar).
  ALE::Obj<real_section_type> _parameters;

  /// Spatial database for final slip.
  spatialdata::spatialdb::SpatialDB* _dbFinalSlip;

  /// Spatial database for slip time.
  spatialdata::spatialdb::SpatialDB* _dbSlipTime;

   /// Spatial database for rise time.
  spatialdata::spatialdb::SpatialDB* _dbRiseTime;

  int _spaceDim; ///< Spatial dimension for slip field.

}; // class LiuCosSlipFn

#include "LiuCosSlipFn.icc" // inline methods

#endif // pylith_faults_liucosslipfn_hh


// End of file 
