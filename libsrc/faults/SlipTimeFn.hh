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

  namespace feassemble {
    class ParameterManager; // HOLDSA ParameterManager
  } // feassemble
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

  /** Create copy of fault.
   *
   * @returns Copy of fault.
   */
  virtual
  SlipTimeFn* clone(void) const = 0;

  /** Initialize slip time function.
   *
   * @param dbSlip Spatial database for slip.
   * @param dbSlipTime Spatial database for slip initiation time.
   * @param dbPeakRate Spatial database for peak slip rate.
   */
  virtual
  void initialize(const ALE::Obj<Mesh>& mesh,
		  const spatialdata::geocoords::CoordSys* cs,
		  const std::set<Mesh::point_type>& vertices) = 0;

  /** Get slip on fault surface at time t.
   *
   * @param t Time t.
   * @param vertices Vertices on fault surface.
   */
  virtual
  const ALE::Obj<real_section_type>& slip(const double t,
			      const std::set<Mesh::point_type>& vertices) = 0;

  // PROTECTED METHODS //////////////////////////////////////////////////
protected :

  /** Copy constructor.
   *
   * @param f SlipTimeFn to copy
   */
  SlipTimeFn(const SlipTimeFn& f);

  // NOT IMPLEMENTED ////////////////////////////////////////////////////
private :

  /// Not implemented
  const SlipTimeFn& operator=(const SlipTimeFn& f);

  // PROTECTED MEMBERS //////////////////////////////////////////////////
protected :

  /// Parameters for slip time function
  feassemble::ParameterManager* _parameters;

}; // class SlipTimeFn

#endif // pylith_faults_sliptimefn_hh


// End of file 
