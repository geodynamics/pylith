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
 */

#if !defined(pylith_faults_faultcohesivekin_hh)
#define pylith_faults_faultcohesivekin_hh

// Include directives ---------------------------------------------------
#include "FaultCohesive.hh" // ISA FaultCohesive

#include "pylith/topology/SubMesh.hh" // ISA Integrator<Quadrature<SubMesh> >
#include "pylith/feassemble/Quadrature.hh" // ISA Integrator<Quadrature>
#include "pylith/feassemble/Integrator.hh" // ISA Integrator

#include <string> // HASA std::string

// FaultCohesiveKin -----------------------------------------------------
/**
 * @brief C++ implementation for a fault surface with kinematic
 * (prescribed) slip implemented with cohesive elements.
 *
 * Fault boundary condition is specified using Lagrange
 * multipliers. The constraints are associated with "constraint"
 * vertices which sit between the pair of vertices on each side of the
 * fault. 
 *
 * The ordering of vertices in a cohesive cell is the vertices on the
 * "negative" side of the fault, the corresponding entries on the
 * "positive" side of the fault, and then the corresponding constraint
 * vertices.
 *
 * The system without Lagrange multipliers is
 *
 * [A(t+dt)]{u(t+dt)} = {b(t+dt)}
 *
 * With Lagrange multipliers this system becomes
 *
 * [A(t+dt) C^T ]{ u(t+dt) } = {b(t+dt)}
 * [ C      0   ]{ L(t+dt) }   {D(t+dt)}
 *
 * where C is the matrix of Lagrange constraints, L is the vector of
 * Lagrange multiplies (internal forces in this case), and D is the
 * fault slip.
 *
 * We solve for the increment in the displacement field, so we rewrite
 * the system as
 *
 * [A(t+dt) C^T ]{ du(t) } = {b(t+dt)} - [A(t+dt) C^T ]{ u(t) }
 * [ C      0   ]{ dL(t) }   {D(t+dt)}   [ C      0   ]{ L(t) }
 * 
 * We form the residual as
 *
 * {r(t+dt)} = {b(t+dt)} - [A(t+dt) C^T ]{ u(t)+du(t) }
 *             {D(t+dt)}   [ C      0   ]{ L(t)+dL(t) }
 * 
 * The terms in the residual contributing to the DOF at the Lagrange
 * vertices are
 *
 * {r(t+dt)} = {D(t+dt)} - [C]{u(t)+dt(t)}
 *
 * The first term, {D(t+dt)}, does not involve integration over the
 * cohesive cells, so it does not require assembling over cohesive
 * cells or processors. We compute the term in
 * integrateResidualAssembled().
 *
 * The term in the residual contributing to the DOF at the
 * non-Lagrange vertices of the cohesive cells is
 *
 * {r(t+dt)} = -[C]^T{L(t)+dL(t)}
 *
 * We integrate the Lagrange multiplier term and the relative
 * displacement term over the cohesive cells, because this introduces
 * weighting of the orientation of the fault for the direction of slip
 * at the vertices of the cohesive cells.
 */
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

  /// Deallocate PETSc and local data structures.
  virtual
  void deallocate(void);
  
  /** Set kinematic earthquake sources.
   *
   * @param names Array of kinematic earthquake source names.
   * @param numNames Number of earthquake sources.
   * @param sources Array of kinematic earthquake sources.
   * @param numSources Number of earthquake sources.
   */
  void eqsrcs(const char* const* names,
	      const int numNames,
	      EqKinSrc** sources,
	      const int numSources);

  /** Initialize fault. Determine orientation and setup boundary
   * condition parameters.
   *
   * @param mesh Finite-element mesh.
   * @param upDir Direction perpendicular to along-strike direction that is 
   *   not collinear with fault normal (usually "up" direction but could 
   *   be up-dip direction; only applies to fault surfaces in a 3-D domain).
   * @param normalDir General preferred direction for fault normal
   *   (used to pick which of two possible normal directions for
   *   interface; only applies to fault surfaces in a 3-D domain).
   */
  void initialize(const topology::Mesh& mesh,
		  const double upDir[3],
		  const double normalDir[3]);

  /** Split solution field for separate preconditioning.
   *
   * @param field Solution field.
   */
  void splitField(topology::Field<topology::Mesh>* field);

  /** Integrate contributions to residual term (r) for operator that
   * require assembly across processors.
   *
   * @param residual Field containing values for residual
   * @param t Current time
   * @param fields Solution fields
   */
  void integrateResidual(const topology::Field<topology::Mesh>& residual,
			 const double t,
			 topology::SolutionFields* const fields);

  /** Integrate contributions to residual term (r) for operator that
   * do not require assembly across cells, vertices, or processors.
   *
   * @param residual Field containing values for residual
   * @param t Current time
   * @param fields Solution fields
   */
  void integrateResidualAssembled(const topology::Field<topology::Mesh>& residual,
				  const double t,
				  topology::SolutionFields* const fields);

  /** Integrate contributions to Jacobian matrix (A) associated with
   * operator that do not require assembly across cells, vertices, or
   * processors.
   *
   * @param jacobian Sparse matrix
   * @param t Current time
   * @param fields Solution fields
   */
  void integrateJacobianAssembled(topology::Jacobian* jacobian,
				  const double t,
				  topology::SolutionFields* const fields);

  /** Integrate contributions to Jacobian matrix (A) associated with
   * operator that do not require assembly across cells, vertices, or
   * processors.
   *
   * @param jacobian Diagonal Jacobian matrix as a field.
   * @param t Current time
   * @param fields Solution fields
   */
  void integrateJacobianAssembled(topology::Field<topology::Mesh>* jacobian,
				  const double t,
				  topology::SolutionFields* const fields);

  /** Update state variables as needed.
   *
   * @param t Current time
   * @param fields Solution fields
   * @param mesh Finite-element mesh
   */
  void updateStateVars(const double t,
		       topology::SolutionFields* const fields);

  /** Adjust solution from solver with lumped Jacobian to match Lagrange
   *  multiplier constraints.
   *
   * @param solution Solution field.
   * @param jacobian Jacobian of the system.
   * @param residual Residual field.
   */
  void adjustSolnLumped(topology::SolutionFields* fields,
			const topology::Field<topology::Mesh>& jacobian);

  /** Verify configuration is acceptable.
   *
   * @param mesh Finite-element mesh
   */
  void verifyConfiguration(const topology::Mesh& mesh) const;

  /** Get vertex field associated with integrator.
   *
   * @param name Name of cell field.
   * @param fields Solution fields.
   * @returns Vertex field.
   */
  const topology::Field<topology::SubMesh>&
  vertexField(const char* name,
	      const topology::SolutionFields* fields =0);

  /** Get cell field associated with integrator.
   *
   * @param name Name of cell field.
   * @param fields Solution fields.
   * @returns Cell field.
   */
  const topology::Field<topology::SubMesh>&
  cellField(const char* name,
	    const topology::SolutionFields* fields =0);

  /** Cohesive cells use Lagrange multiplier constraints?
   *
   * @returns True if implementation using Lagrange multiplier
   * constraints, false otherwise.
   */
  bool useLagrangeConstraints(void) const;

  /** Get fields associated with fault.
   *
   * @returns Fields associated with fault.
   */
  const topology::Fields<topology::Field<topology::SubMesh> >*
  fields(void) const;

  // PRIVATE METHODS ////////////////////////////////////////////////////
private :

  /// Initialize logger.
  void _initializeLogger(void);
  
  /** Calculate orientation at fault vertices.
   *
   * @param upDir Direction perpendicular to along-strike direction that is 
   *   not collinear with fault normal (usually "up" direction but could 
   *   be up-dip direction; only applies to fault surfaces in a 3-D domain).
   * @param normalDir General preferred direction for fault normal
   *   (used to pick which of two possible normal directions for
   *   interface; only applies to fault surfaces in a 3-D domain).
   */
  void _calcOrientation(const double upDir[3],
			const double normalDir[3]);

  /// Calculate fault area field.
  void _calcArea(void);

  /** Compute change in tractions on fault surface using solution.
   *
   * @param tractions Field for tractions.
   * @param mesh Finite-element mesh for domain
   * @param solution Solution over domain
   */
  void _calcTractionsChange(topology::Field<topology::SubMesh>* tractions,
			    const topology::Field<topology::Mesh>& solution);

  /// Allocate buffer for vector field.
  void _allocateBufferVectorField(void);

  /// Allocate buffer for scalar field.
  void _allocateBufferScalarField(void);

  // PRIVATE TYPEDEFS ///////////////////////////////////////////////////
private :

  typedef std::map<std::string, EqKinSrc*> srcs_type;

  // PRIVATE MEMBERS ////////////////////////////////////////////////////
private :

  srcs_type _eqSrcs; ///< Array of kinematic earthquake sources.

  // NOT IMPLEMENTED ////////////////////////////////////////////////////
private :

  /// Not implemented
  FaultCohesiveKin(const FaultCohesiveKin&);

  /// Not implemented
  const FaultCohesiveKin& operator=(const FaultCohesiveKin&);

}; // class FaultCohesiveKin

#include "FaultCohesiveKin.icc" // inline methods

#endif // pylith_faults_faultcohesivekin_hh


// End of file 
