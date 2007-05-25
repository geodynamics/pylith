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
 * @brief C++ implementation for a fault surface with kinematic
 * (prescribed) slip implemented with cohesive elements.
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

#include "FaultCohesive.hh" // ISA FaultCohesive
#include "pylith/feassemble/Constraint.hh" // ISA Constraint
#include "pylith/feassemble/Integrator.hh" // ISA Integrator

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
class pylith::faults::FaultCohesiveKin : public FaultCohesive,
					 public feassemble::Integrator,
					 public feassemble::Constraint
{ // class FaultCohesiveKin
  friend class TestFaultCohesiveKin; // unit testing

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Default constructor.
  FaultCohesiveKin(void);

  /// Destructor.
  virtual
  ~FaultCohesiveKin(void);

  /** Set kinematic earthquake source.
   *
   * @param src Kinematic earthquake source.
   */
  void eqsrc(EqKinSrc* src);

  /** Initialize fault. Determine orientation and setup boundary
   * condition parameters.
   *
   * @param mesh PETSc mesh
   * @param cs Coordinate system for mesh
   * @param upDir Direction perpendicular to along-strike direction that is 
   *   not collinear with fault normal (usually "up" direction but could 
   *   be up-dip direction).
   */
  void initialize(const ALE::Obj<ALE::Mesh>& mesh,
		  const spatialdata::geocoords::CoordSys* cs,
		  const double_array& upDir);

  /** Integrate contributions to residual term (r) for operator.
   *
   * @param residual Field containing values for residual
   * @param fields Solution fields
   * @param mesh Finite-element mesh
   */
  void integrateResidual(const ALE::Obj<real_section_type>& residual,
			 topology::FieldsManager* const fields,
			 const ALE::Obj<Mesh>& mesh);

  /** Integrate contributions to Jacobian matrix (A) associated with
   * operator.
   *
   * @param mat Sparse matrix
   * @param fields Solution fields
   * @param mesh Finite-element mesh
   */
  void integrateJacobian(PetscMat* mat,
			 topology::FieldsManager* const fields,
			 const ALE::Obj<Mesh>& mesh);

  /** Set number of degrees of freedom that are constrained at points
   * in field.
   *
   * @param field Solution field
   * @param mesh PETSc mesh
   */
  void setConstraintSizes(const ALE::Obj<real_section_type>& field,
			  const ALE::Obj<ALE::Mesh>& mesh);

  /** Set which degrees of freedom are constrained at points in field.
   *
   * @param field Solution field
   * @param mesh PETSc mesh
   */
  void setConstraints(const ALE::Obj<real_section_type>& field,
		      const ALE::Obj<ALE::Mesh>& mesh);

  /** Set field.
   *
   * @param t Current time.
   * @param disp Displacement field at time t.
   * @param mesh Finite-element mesh.
   */
  void setField(const double t,
		const ALE::Obj<real_section_type>& disp,
		const ALE::Obj<Mesh>& mesh);
  
  // PROTECTED METHODS //////////////////////////////////////////////////
protected :

  /** Cohesive cells use Lagrange multiplier constraints?
   *
   * @returns True if implementation using Lagrange multiplier
   * constraints, false otherwise.
   */
  bool _useLagrangeConstraints(void) const;

  // NOT IMPLEMENTED ////////////////////////////////////////////////////
private :

  /// Not implemented
  FaultCohesiveKin(const FaultCohesiveKin& m);

  /// Not implemented
  const FaultCohesiveKin& operator=(const FaultCohesiveKin& m);

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
