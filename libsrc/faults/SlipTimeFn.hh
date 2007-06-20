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

  /** Get slip increment on fault surface between time t0 and t1.
   *
   * @param t0 Time t.
   * @param t1 Time t+dt.
   * @param vertices Vertices where slip will be prescribed.
   */
  virtual
  const ALE::Obj<real_section_type>& slipIncr(const double t0,
					      const double t1,
					      const std::set<Mesh::point_type>& vertices) = 0;

  // NOT IMPLEMENTED ////////////////////////////////////////////////////
private :

  /// Not implemented.
  SlipTimeFn(const SlipTimeFn& f);

  /// Not implemented
  const SlipTimeFn& operator=(const SlipTimeFn& f);

  // PROTECTED MEMBERS //////////////////////////////////////////////////
protected :

  /// Parameters for slip time function
  topology::FieldsManager* _parameters;

}; // class SlipTimeFn

#endif // pylith_faults_sliptimefn_hh


// End of file 
