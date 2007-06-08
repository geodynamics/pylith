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

/** @file libsrc/faults/FaultCohesive.hh
 *
 * @brief C++ abstract base class for a fault surface implemented with
 * cohesive elements.
 */

#if !defined(pylith_faults_faultcohesive_hh)
#define pylith_faults_faultcohesive_hh

#include "Fault.hh" // ISA Fault

#include "pylith/utils/sievefwd.hh" // HOLDSA PETSc Mesh
#include "pylith/utils/petscfwd.h" // USES PetscMat

/// Namespace for pylith package
namespace pylith {
  namespace faults {
    class FaultCohesive;
    class TestFaultCohesive; // unit testing
  } // faults
} // pylith

/// C++ abstract base class for a fault surface implemented with
/// cohesive elements.
class pylith::faults::FaultCohesive : public Fault
{ // class FaultCohesive
  friend class TestFaultCohesive; // unit testing

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Default constructor.
  FaultCohesive(void);

  /// Destructor.
  virtual
  ~FaultCohesive(void);

  /** Adjust mesh topology for fault implementation.
   *
   * @param mesh PETSc mesh
   */
  void adjustTopology(const ALE::Obj<ALE::Mesh>& mesh);

  // PROTECTED TYPEDEFS /////////////////////////////////////////////////
protected :

  /// Function type for orientation methods.
  typedef void (*orient_fn_type)(double_array*, 
				 const double_array&,
				 const double,
				 const double_array&);

  // PROTECTED METHODS //////////////////////////////////////////////////
protected :

  /** Cohesive cells use Lagrange multiplier constraints?
   *
   * @returns True if implementation using Lagrange multiplier
   * constraints, false otherwise.
   */
  virtual
  bool _useLagrangeConstraints(void) const = 0;

  /** Compute weighted orientation of fault for cohesive cell between
   * 1-D elements. Orientation is either at vertices or quadrature
   * points, depending on whether the arguments have been evaluated at
   * the vertices or quadrature points.
   *
   * The orientation is returned as an array of direction cosines.
   *
   * size = spaceDim*spaceDim
   * index = iDir*spaceDim + iComponent
   *
   * @param orientation Array of direction cosines.
   * @param jacobian Jacobian matrix at point.
   * @param jacobianDet Determinant of Jacobian matrix at point.
   * @param upDir Direction perpendicular to along-strike direction that is 
   *   not collinear with fault normal (usually "up" direction but could 
   *   be up-dip direction).
   */
  static
  void _orient1D(double_array* orientation,
		 const double_array& jacobian,
		 const double jacobianDet,
		 const double_array& upDir);
		
  /** Compute weighted orientation of fault for cohesive cell between
   * 2-D elements. Orientation is either at vertices or quadrature
   * points, depending on whether the arguments have been evaluated at
   * the vertices or quadrature points.
   *
   * The orientation is returned as an array of direction cosines.
   *
   * size = spaceDim*spaceDim
   * index = iDir*spaceDim + iComponent
   *
   * @param orientation Array of direction cosines.
   * @param jacobian Jacobian matrix at point.
   * @param jacobianDet Determinant of Jacobian matrix at point.
   * @param upDir Direction perpendicular to along-strike direction that is 
   *   not collinear with fault normal (usually "up" direction but could 
   *   be up-dip direction).
   */
  static 
  void _orient2D(double_array* orientation,
		 const double_array& jacobian,
		 const double jacobianDet,
		 const double_array& upDir);
		
  /** Compute weighted orientation of fault for cohesive cell between
   * 3-D elements. Orientation is either at vertices or quadrature
   * points, depending on whether the arguments have been evaluated at
   * the vertices or quadrature points.
   *
   * The orientation is returned as an array of direction cosines.
   *
   * size = spaceDim*spaceDim
   * index = iDir*spaceDim + iComponent
   *
   * @param orientation Array of direction cosines.
   * @param jacobian Jacobian matrix at point.
   * @param jacobianDet Determinant of Jacobian matrix at point.
   * @param upDir Direction perpendicular to along-strike direction that is 
   *   not collinear with fault normal (usually "up" direction but could 
   *   be up-dip direction).
   */
  static
  void _orient3D(double_array* orientation,
		 const double_array& jacobian,
		 const double jacobianDet,
		 const double_array& upDir);
		
  // NOT IMPLEMENTED ////////////////////////////////////////////////////
private :

  /// Not implemented
  FaultCohesive(const FaultCohesive& m);

  /// Not implemented
  const FaultCohesive& operator=(const FaultCohesive& m);

// PROTECTED MEMBERS ////////////////////////////////////////////////////
protected :

  ALE::Obj<ALE::Mesh>* _faultMesh;

}; // class FaultCohesive

#endif // pylith_faults_faultcohesive_hh


// End of file 
