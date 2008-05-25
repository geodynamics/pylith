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

/** @file libsrc/faults/EqKinSrc.hh
 *
 * @brief C++ object for managing parameters for a kinematic
 * earthquake source.
 *
 * EqKinSrc is responsible for providing the value of slip at time t
 * over a fault surface.
 */

#if !defined(pylith_faults_eqkinsrc_hh)
#define pylith_faults_eqkinsrc_hh

#include "pylith/utils/sievetypes.hh" // USES PETSc Mesh

/// Namespace for pylith package
namespace pylith {
  namespace faults {
    class EqKinSrc;
    class TestEqKinSrc; // unit testing

    class SlipTimeFn; // HOLDSA SlipTimeFn
  } // faults
} // pylith

/// Namespace for spatialdata package
namespace spatialdata {
  namespace spatialdb {
    class SpatialDB;
  } // spatialdb

  namespace geocoords {
    class CoordSys;
  } // geocoords
} // spatialdata

/// C++ oject for managing parameters for a kinematic earthquake source.
class pylith::faults::EqKinSrc
{ // class EqKinSrc
  friend class TestEqKinSrc; // unit testing

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Default constructor.
  EqKinSrc(void);

  /// Destructor.
  virtual
  ~EqKinSrc(void);

  /** Set origin time for earthquake source.
   *
   * @param value Origin time for earthquake source.
   */
  void originTime(const double value);

  /** Get origin time for earthquake source.
   *
   * @returns Origin time for earthquake source.
   */
  double originTime(void) const;

  /** Set slip time function.
   *
   * @param slipfn Slip time function.
   */
  void slipfn(SlipTimeFn* slipfn);

  /** Initialize slip time function.
   *
   * @param faultMesh Finite-element mesh of fault.
   * @param cs Coordinate system for mesh
   */
  virtual
  void initialize(const ALE::Obj<Mesh>& faultMesh,
		  const spatialdata::geocoords::CoordSys* cs);

  /** Get slip on fault surface at time t.
   *
   * @param t Time t.
   * @param faultMesh Finite-element mesh of fault.
   */
  virtual
  const ALE::Obj<real_section_type>& slip(const double t,
					  const ALE::Obj<Mesh>& faultMesh);

  /** Get increment of slip on fault surface between time t0 and t1.
   *
   * @param t Time t.
   * @param faultMesh Finite-element mesh of fault.
   */
  virtual
  const ALE::Obj<real_section_type>& slipIncr(
			      const double t0,
			      const double t1,
			      const ALE::Obj<Mesh>& faultMesh);

  /** Get final slip.
   *
   * @returns Final slip.
   */
  virtual
  ALE::Obj<real_section_type> finalSlip(void);

  /** Get time when slip begins at each point.
   *
   * @returns Time when slip begins.
   */
  virtual
  ALE::Obj<real_section_type> slipTime(void);

  // NOT IMPLEMENTED ////////////////////////////////////////////////////
private :

  /// Not implemented
  EqKinSrc(const EqKinSrc& s);

  /// Not implemented
  const EqKinSrc& operator=(const EqKinSrc& s);

  // PRIVATE MEMBERS ////////////////////////////////////////////////////
private :

  double _originTime; ///< Origin time for earthquake source
  SlipTimeFn* _slipfn; ///< Slip time function

}; // class EqKinSrc

#endif // pylith_faults_eqkinsrc_hh


// End of file 
