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

// Include directives ---------------------------------------------------
#include "faultsfwd.hh" // forward declarations

#include "pylith/topology/topologyfwd.hh" // USES Field<SubMesh>

#include "spatialdata/units/unitsfwd.hh" // USES Nondimensional

// EqKinSrc -------------------------------------------------------------
class pylith::faults::EqKinSrc
{ // class EqKinSrc
  friend class TestEqKinSrc; // unit testing

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Default constructor.
  EqKinSrc(void);

  /// Destructor.
  ~EqKinSrc(void);

  /// Deallocate PETSc and local data structures.
  virtual
  void deallocate(void);
  
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
   * @param normalizer Nondimensionalization of scales.
   */
  void initialize(const topology::SubMesh& faultMesh,
		  const spatialdata::units::Nondimensional& normalizer);

  /** Get slip on fault surface at time t.
   *
   * @param slipField Slip field over fault mesh.
   * @param t Time t.
   */
  void slip(topology::Field<topology::SubMesh>* const slipField,
	    const double t);

  /** Get increment of slip on fault surface between time t0 and t1.
   *
   * @param slipField Slip increment field over fault mesh.
   * @param t0 Time for start of slip increment.
   * @param t1 Time for end of slip increment.
   */
  void slipIncr(topology::Field<topology::SubMesh>* const slipField,
		const double t0,
		const double t1);

  /** Get final slip.
   *
   * @returns Final slip.
   */
  const topology::Field<topology::SubMesh>& finalSlip(void) const;

  /** Get time when slip begins at each point.
   *
   * @returns Time when slip begins.
   */
  const topology::Field<topology::SubMesh>& slipTime(void) const;

  // NOT IMPLEMENTED ////////////////////////////////////////////////////
private :

  EqKinSrc(const EqKinSrc&); ///< Not implemented
  const EqKinSrc& operator=(const EqKinSrc&); ///< Not implemented

  // PRIVATE MEMBERS ////////////////////////////////////////////////////
private :

  double _originTime; ///< Origin time for earthquake source
  SlipTimeFn* _slipfn; ///< Slip time function

}; // class EqKinSrc

#endif // pylith_faults_eqkinsrc_hh


// End of file 
