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

  namespace topology {
    class FieldsManager; // HOLDSA FieldsManager
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
   * @param mesh Finite-element mesh.
   * @param faultMesh Finite-element mesh of fault.
   * @param vertices Vertices where slip will be prescribed.
   * @param cs Coordinate system for mesh
   */
  virtual
  void initialize(const ALE::Obj<Mesh>& mesh,
		  const ALE::Obj<Mesh>& faultMesh,
		  const std::set<Mesh::point_type>& vertices,
		  const spatialdata::geocoords::CoordSys* cs) = 0;

  /** Get slip on fault surface at time t.
   *
   * @param t Time t.
   * @param vertices Vertices where slip will be prescribed.
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
  topology::FieldsManager* _parameters;

}; // class SlipTimeFn

#endif // pylith_faults_sliptimefn_hh


// End of file 
