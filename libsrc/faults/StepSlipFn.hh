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

/** @file libsrc/faults/StepSlipFn.hh
 *
 * @brief C++ implementation of a step-function slip time function.
 *
 * Slip time function is a step function with slip beginning at time t0.
 *
 * Normalized slip = 1 if t >= t0, 0 otherwise
 */

#if !defined(pylith_faults_stepslipfn_hh)
#define pylith_faults_stepslipfn_hh

#include "SlipTimeFn.hh"

/// Namespace for pylith package
namespace pylith {
  namespace faults {
    class StepSlipFn;
    class TestStepSlipFn; // unit testing
  } // faults
} // pylith

/// Namespace for spatialdata package
namespace spatialdata {
  namespace spatialdb {
    class SpatialDB;
  } // spatialdb
} // spatialdata

/// C++ implementation of Step slip time function.
class pylith::faults::StepSlipFn : public SlipTimeFn
{ // class StepSlipFn
  friend class TestStepSlipFn; // unit testing

// PUBLIC METHODS ///////////////////////////////////////////////////////
public :

  /// Default constructor.
  StepSlipFn(void);

  /// Destructor.
  ~StepSlipFn(void);

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
  StepSlipFn(const StepSlipFn& m);

  /// Not implemented
  const StepSlipFn& operator=(const StepSlipFn& f);

// PRIVATE MEMBERS //////////////////////////////////////////////////////
private :

  /// Parameters for Step slip time function.
  /// Final slip (vector), slip time (scalar).
  ALE::Obj<real_section_type> _parameters;

  /// Spatial database for final slip
  spatialdata::spatialdb::SpatialDB* _dbFinalSlip;

  /// Spatial database for slip time
  spatialdata::spatialdb::SpatialDB* _dbSlipTime;

  int _spaceDim; ///< Spatial dimension for slip field.

}; // class StepSlipFn

#include "StepSlipFn.icc" // inline methods

#endif // pylith_faults_stepslipfn_hh


// End of file 
