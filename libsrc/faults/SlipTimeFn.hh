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

/** @file libsrc/faults/SlipTimeFn.hh
 *
 * @brief C++ abstract base class for kinematic slip time function.
 *
 * Interface definition for slip time function.
 */

#if !defined(pylith_faults_sliptimefn_hh)
#define pylith_faults_sliptimefn_hh

// Include directives ---------------------------------------------------
#include "faultsfwd.hh" // forward declarations

#include "pylith/topology/topologyfwd.hh" // USES Field<SubMesh>

#include "spatialdata/units/unitsfwd.hh" // USES Nondimensional
#include "spatialdata/spatialdb/spatialdbfwd.hh" // USES SpatialDB

// SlipTimeFn -----------------------------------------------------------
class pylith::faults::SlipTimeFn
{ // class SlipTimeFn
  friend class TestSlipTimeFn; // unit testing

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Default constructor.
  SlipTimeFn(void);

  /// Destructor.
  virtual
  ~SlipTimeFn(void);

  /** Initialize slip time function.
   *
   * @param faultMesh Finite-element mesh of fault.
   * @param cs Coordinate system for mesh
   * @param normalizer Nondimensionalization of scales.
   * @param originTime Origin time for earthquake source.
   */
  virtual
  void initialize(const topology::SubMesh& faultMesh,
		  const spatialdata::units::Nondimensional& normalizer,
		  const double originTime =0.0) = 0;

  /** Get slip on fault surface at time t.
   *
   * @param slipField Slip field over fault surface.
   * @param t Time t.
   *
   * @returns Slip vector as left-lateral/reverse/normal.
   */
  virtual
  void slip(topology::Field<topology::SubMesh>* const slipField,
	    const double t) = 0;
  
  /** Get slip increment on fault surface between time t0 and t1.
   *
   * @param slipField Slip field over fault surface.
   * @param t0 Time t.
   * @param t1 Time t+dt.
   * 
   * @returns Increment in slip vector as left-lateral/reverse/normal.
   */
  virtual
  void slipIncr(topology::Field<topology::SubMesh>* slipField,
		const double t0,
		const double t1) = 0;

  /** Get final slip.
   *
   * @returns Final slip.
   */
  virtual
  const topology::Field<topology::SubMesh>& finalSlip(void) = 0;

  /** Get time when slip begins at each point.
   *
   * @returns Time when slip begins.
   */
  virtual
  const topology::Field<topology::SubMesh>& slipTime(void) = 0;

  // NOT IMPLEMENTED ////////////////////////////////////////////////////
private :

  SlipTimeFn(const SlipTimeFn&); ///< Not implemented
  const SlipTimeFn& operator=(const SlipTimeFn&); ///< Not implemented

}; // class SlipTimeFn

#endif // pylith_faults_sliptimefn_hh


// End of file 
