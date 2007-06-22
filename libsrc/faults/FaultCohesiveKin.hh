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
 * one side of the fault, the corresponding entries on the other side
 * of the fault, and then the corresponding constraint vertices.
 *
 * [ K  aC^T ] [ U  ] = [ Fe ]
 * [ C   0   ] [ Fi ] = [ D  ]
 *
 * where K is the stiffness matrix, C is the matrix of Lagrange
 * constraints, U is the displacement field, Fe is the vector of
 * external forces, Fi is the vector of Lagrange multipers (forces), D
 * is the fault slip, and "a" is the conditioning value (taken to be
 * the shear modulus).
 */

#if !defined(pylith_faults_faultcohesivekin_hh)
#define pylith_faults_faultcohesivekin_hh

#include "FaultCohesive.hh" // ISA FaultCohesive
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
					 public feassemble::Integrator
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
   *   be up-dip direction; only applies to fault surfaces in a 3-D domain).
   * @param normalDir General preferred direction for fault normal
   *   (used to pick which of two possible normal directions for
   *   interface; only applies to fault surfaces in a 3-D domain).
   * @param matDB Database of bulk elastic properties for fault region
   *   (used to improve conditioning of Jacobian matrix)
   */
  void initialize(const ALE::Obj<ALE::Mesh>& mesh,
		  const spatialdata::geocoords::CoordSys* cs,
		  const double_array& upDir,
		  const double_array& normalDir,
		  spatialdata::spatialdb::SpatialDB* matDB);

  /** Integrate contributions to residual term (r) for operator.
   *
   * @param residual Field containing values for residual
   * @param t Current time
   * @param fields Solution fields
   * @param mesh Finite-element mesh
   */
  void integrateResidual(const ALE::Obj<real_section_type>& residual,
			 const double t,
			 topology::FieldsManager* const fields,
			 const ALE::Obj<Mesh>& mesh);

  /** Integrate contributions to Jacobian matrix (A) associated with
   * operator.
   *
   * @param mat Sparse matrix
   * @param t Current time
   * @param fields Solution fields
   * @param mesh Finite-element mesh
   */
  void integrateJacobian(PetscMat* mat,
			 const double t,
			 topology::FieldsManager* const fields,
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

  /// Pseudo-stiffness for scaling constraint information to improve
  /// conditioning of Jacobian matrix
  ALE::Obj<real_section_type> _pseudoStiffness;

  /// Orientation of fault surface at vertices (fiber dimension is
  /// nonzero only at constraint vertices)
  ALE::Obj<real_section_type> _orientation;

  /// Fault vertices associated with constraints
  std::set<Mesh::point_type> _constraintVert;

  /// Label of cell used to compute Jacobian for each constraint vertex (must
  /// prevent overlap so that only allow 1 cell will contribute for
  /// each vertex).
  ALE::Obj<int_section_type> _constraintCell;



}; // class FaultCohesiveKin

#include "FaultCohesiveKin.icc" // inline methods

#endif // pylith_faults_faultcohesivekin_hh


// End of file 
