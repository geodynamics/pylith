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

#include "pylith/utils/sievetypes.hh" // USES PETSc Mesh

/// Namespace for pylith package
namespace pylith {
  namespace faults {
    class SlipTimeFn;
    class TestSlipTimeFn; // unit testing
  } // faults
} // pylith

/// Namespace for spatialdata package
namespace spatialdata {
  namespace geocoords {
    class CoordSys;
  } // geocoords
} // spatialdata

/// C++ abstract base class for Fault object.
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
   * @param originTime Origin time for earthquake source.
   */
  virtual
  void initialize(const ALE::Obj<Mesh>& faultMesh,
		  const spatialdata::geocoords::CoordSys* cs,
		  const double originTime =0.0) = 0;

  /** Get slip on fault surface at time t.
   *
   * @param t Time t.
   * @param faultMesh Mesh over fault surface.
   *
   * @returns Slip vector as left-lateral/reverse/normal.
   */
  virtual
  const ALE::Obj<real_section_type>&
  slip(const double t,
       const ALE::Obj<Mesh>& faultMesh) = 0;
  
  /** Get slip increment on fault surface between time t0 and t1.
   *
   * @param t0 Time t.
   * @param t1 Time t+dt.
   * @param faultMesh Mesh over fault surface.
   * 
   * @returns Increment in slip vector as left-lateral/reverse/normal.
   */
  virtual
  const ALE::Obj<real_section_type>&
  slipIncr(const double t0,
	   const double t1,
	   const ALE::Obj<Mesh>& faultMesh) = 0;

  /** Get final slip.
   *
   * @returns Final slip.
   */
  virtual
  ALE::Obj<real_section_type> finalSlip(void) = 0;

  /** Get time when slip begins at each point.
   *
   * @returns Time when slip begins.
   */
  virtual
  ALE::Obj<real_section_type> slipTime(void) = 0;

  // NOT IMPLEMENTED ////////////////////////////////////////////////////
private :

  /// Not implemented.
  SlipTimeFn(const SlipTimeFn& f);

  /// Not implemented
  const SlipTimeFn& operator=(const SlipTimeFn& f);

}; // class SlipTimeFn

#endif // pylith_faults_sliptimefn_hh


// End of file 
