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

/** @file libsrc/faults/FaultCohesiveKin.hh
 *
 * @brief C++ object implementing a fault surface with a kinematic
 * earthquake source.
 *
 * Fault boundary condition is specified using Lagrange
 * multipliers. The constraints are associated with "constraint"
 * vertices which sit between the pair of vertices on each side of the
 * fault. 
 *
 * The ordering of vertices in a cohesive cell is the vertices on the
 * POSITIVE/NEGATIVE (CHECK WHICH IT IS) side of the fault, the
 * corresponding entries on the other side of the fault, and then the
 * corresponding constraint vertices.
 */

#if !defined(pylith_faults_faultcohesivekin_hh)
#define pylith_faults_faultcohesivekin_hh

#include "FaultCohesive.hh"

/// Namespace for pylith package
namespace pylith {
  namespace faults {
    class FaultCohesiveKin;
    class TestFaultCohesiveKin; // unit testing

    class EqKinSrc; // HOLDSA EqKinSrc
  } // faults
} // pylith

/// C++ implementation for a fault surface with kinematic (prescribed)
/// slip implemented with cohesive elements.
class pylith::faults::FaultCohesiveKin : public FaultCohesive
{ // class FaultCohesiveKin
  friend class TestFaultCohesiveKin; // unit testing

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Default constructor.
  FaultCohesiveKin(void);

  /// Destructor.
  virtual
  ~FaultCohesiveKin(void);

  /** Create copy of fault.
   *
   * @returns Copy of fault.
   */
  Fault* clone(void) const;  

  /** Set kinematic earthquake source.
   *
   * @param src Kinematic earthquake source.
   */
  void eqsrc(EqKinSrc* src);

  /** Initialize fault. Determine orientation and setup boundary
   * condition parameters.
   *
   * @param mesh PETSc mesh
   * @param upDir Direction perpendicular to along-strike direction that is 
   *   not collinear with fault normal (usually "up" direction but could 
   *   be up-dip direction).
   */
  void initialize(const ALE::Obj<ALE::Mesh>& mesh,
		  const double_array& upDir);

  /** Integrate contribution of cohesive cells to residual term.
   *
   * @param residual Residual field (output)
   * @param disp Displacement field at time t
   * @param mesh Finite-element mesh
   */
  void integrateResidual(const ALE::Obj<real_section_type>& residual,
			 const ALE::Obj<real_section_type>& disp,
			 const ALE::Obj<Mesh>& mesh);

  /** Compute Jacobian matrix (A) associated with operator.
   *
   * @param mat Sparse matrix
   * @param disp Displacement field
   * @param mesh Finite-element mesh
   */
  void integrateJacobian(PetscMat* mat,
			 const ALE::Obj<real_section_type>& dispT,
			 const ALE::Obj<Mesh>& mesh);
  
  /** Set field.
   *
   * @param t Current time
   * @param disp Displacement field
   * @param mesh Finite-element mesh
   */
  void setField(const double t,
		const ALE::Obj<real_section_type>& disp,
		const ALE::Obj<Mesh>& mesh);
  
  // PROTECTED METHODS //////////////////////////////////////////////////
protected :

  /** Copy constructor.
   *
   * @param m Fault to copy
   */
  FaultCohesiveKin(const FaultCohesiveKin& m);

  /** Cohesive cells use Lagrange multiplier constraints?
   *
   * @returns True if implementation using Lagrange multiplier
   * constraints, false otherwise.
   */
  bool _useLagrangeConstraints(void) const;

  // NOT IMPLEMENTED ////////////////////////////////////////////////////
private :

  /// Not implemented
  const FaultCohesiveKin& operator=(const FaultCohesiveKin& m);

  // PRIVATE METHODS ////////////////////////////////////////////////////
private :

  /** Get size (fiber dimension) of orientation information.
   *
   * @returns Size of orientation information.
   */
  int _orientationSize(void) const;

  // PRIVATE MEMBERS ////////////////////////////////////////////////////
private :

  EqKinSrc* _eqsrc; ///< Kinematic earthquake source information

  /// Orientation of fault surface at vertices (fiber dimension is
  /// nonzero only at constraint vertices)
  ALE::Obj<real_section_type> _orientation;

  /// Fault vertices associated with constraints
  std::set<Mesh::point_type> _constraintVert;

}; // class FaultCohesiveKin

#include "FaultCohesiveKin.icc" // inline methods

#endif // pylith_faults_faultcohesivekin_hh


// End of file 
