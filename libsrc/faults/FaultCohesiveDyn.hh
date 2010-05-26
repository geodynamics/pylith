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

/** @file libsrc/faults/FaultCohesiveDyn.hh
 *
 * @brief C++ implementation for a fault surface with spontaneous
 * (dynamic) slip implemented with cohesive cells.
 */

#if !defined(pylith_faults_faultcohesivedyn_hh)
#define pylith_faults_faultcohesivedyn_hh

// Include directives ---------------------------------------------------
#include "FaultCohesiveLagrange.hh" // ISA FaultCohesiveLagrange

#include "pylith/friction/frictionfwd.hh" // HOLDSA Friction model
#include "pylith/utils/petscfwd.h" // HASA PetscKSP

// FaultCohesiveDyn -----------------------------------------------------
/**
 * @brief C++ implementation for a fault surface with spontaneous
 * (dynamic) slip implemented with cohesive cells.
 *
 */
class pylith::faults::FaultCohesiveDyn : public FaultCohesiveLagrange
{ // class FaultCohesiveDyn
  friend class TestFaultCohesiveDyn; // unit testing

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Default constructor.
  FaultCohesiveDyn(void);

  /// Destructor.
  virtual
  ~FaultCohesiveDyn(void);

  /// Deallocate PETSc and local data structures.
  virtual
  void deallocate(void);
  
  /** Sets the spatial database for the inital tractions.
   *
   * @param db spatial database for initial tractions
   */
  void dbInitialTract(spatialdata::spatialdb::SpatialDB* db);
  
  /** Set the friction (constitutive) model.
   *
   * @param model Fault constutive model.
   */
  void frictionModel(friction::FrictionModel* const model);

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

  /** Integrate contributions to residual term (r) for operator that
   * require assembly processors.
   *
   * Initial tractions (if specified) contribute to the residual like
   * Neumann boundary conditions.
   *
   * @param residual Field containing values for residual
   * @param t Current time
   * @param fields Solution fields
   */
  virtual
  void integrateResidual(const topology::Field<topology::Mesh>& residual,
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

  /** Constrain solution space based on friction.
   *
   * @param fields Solution fields.
   * @param t Current time.
   * @param jacobian Sparse matrix for system Jacobian.
   */
  void constrainSolnSpace(topology::SolutionFields* const fields,
			  const double t,
			  const topology::Jacobian& jacobian);

  /** Adjust solution from solver with lumped Jacobian to match Lagrange
   *  multiplier constraints.
   *
   * @param fields Solution fields.
   * @param jacobian Jacobian of the system.
   */
  void adjustSolnLumped(topology::SolutionFields* fields,
      const topology::Field<topology::Mesh>& jacobian);

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

  // PRIVATE METHODS ////////////////////////////////////////////////////
private :

  /** Get initial tractions using a spatial database.
   */
  void _setupInitialTractions(void);

  /** Compute tractions on fault surface using solution.
   *
   * @param tractions Field for tractions.
   * @param solution Solution over domain
   */
  void _calcTractions(topology::Field<topology::SubMesh>* tractions,
          const topology::Field<topology::Mesh>& solution);

  /** Get initial tractions on fault surface.
   *
   * @param tractions Field for tractions.
   */
  void _getInitialTractions(topology::Field<topology::SubMesh>* tractions);

  /** Update slip rate associated with Lagrange vertex k corresponding
   * to diffential velocity between conventional vertices i and j.
   *
   * @param fields Solution fields.
   */
  void _updateSlipRate(const topology::SolutionFields& fields);

  /** Setup sensitivity problem to compute change in slip given change
   * in Lagrange multipliers.
   *
   * @param jacobian Jacobian matrix for entire domain.
   */
  void _sensitivitySetup(const topology::Jacobian& jacobian);

  /** Update the Jacobian values for the sensitivity solve.
   *
   * @param negativeSide True if solving sensitivity problem for
   * negative side of the fault, false if solving sensitivity problem
   * for positive side of the fault.
   * @param jacobian Jacobian matrix for entire domain.
   * @param fields Solution fields.
   */
  void _sensitivityUpdateJacobian(const bool negativeSide,
                                  const topology::Jacobian& jacobian,
                                  const topology::SolutionFields& fields);

  /** Reform residual for sensitivity problem.
   *
   * @param negativeSide True if solving sensitivity problem for
   * negative side of the fault, false if solving sensitivity problem
   * for positive side of the fault.
   */
  void _sensitivityReformResidual(const bool negativeSide);

  /// Solve sensitivity problem.
  void _sensitivitySolve(void);

  /** Update the solution (displacement increment) values based on
   * the sensitivity solve.
   *
   * @param negativeSide True if solving sensitivity problem for
   * negative side of the fault, false if solving sensitivity problem
   * for positive side of the fault.
   */
  void _sensitivityUpdateSoln(const bool negativeSide);

  /** Solve slip/Lagrange multiplier sensitivity problem for case of lumped Jacobian in 1-D.
   *
   * @param slip Adjustment to slip assoc. w/Lagrange multiplier vertex.
   * @param dLagrangeTpdt Change in Lagrange multiplier.
   * @param jacobianN Jacobian for vertex on - side of the fault.
   * @param jacobianP Jacobian for vertex on + side of the fault.
   */
  void _sensitivitySolveLumped1D(double_array* slip,
                                 const double_array& dLagrangeTpdt,
                                 const double_array& jacobianN,
                                 const double_array& jacobianP);

  /** Solve slip/Lagrange multiplier sensitivity problem for case of lumped Jacobian in 2-D.
   *
   * @param slip Adjustment to slip assoc. w/Lagrange multiplier vertex.
   * @param dLagrangeTpdt Change in Lagrange multiplier.
   * @param jacobianN Jacobian for vertex on - side of the fault.
   * @param jacobianP Jacobian for vertex on + side of the fault.
   */
  void _sensitivitySolveLumped2D(double_array* slip,
                                 const double_array& dLagrangeTpdt,
                                 const double_array& jacobianN,
                                 const double_array& jacobianP);

  /** Solve slip/Lagrange multiplier sensitivity problem for case of lumped Jacobian in 3-D.
   *
   * @param slip Adjustment to slip assoc. w/Lagrange multiplier vertex.
   * @param dLagrangeTpdt Change in Lagrange multiplier.
   * @param jacobianN Jacobian for vertex on - side of the fault.
   * @param jacobianP Jacobian for vertex on + side of the fault.
   */
  void _sensitivitySolveLumped3D(double_array* slip,
                                 const double_array& dLagrangeTpdt,
                                 const double_array& jacobianN,
                                 const double_array& jacobianP);

  /** Constrain solution space with lumped Jacobian in 1-D.
   *
   * @param dLagrangeTpdt Adjustment to Lagrange multiplier.
   * @param slip Slip assoc. w/Lagrange multiplier vertex.
   * @param slipRate Slip rate assoc. w/Lagrange multiplier vertex.
   * @param tractionTpdt Fault traction assoc. w/Lagrange multiplier vertex.
   * @param area Fault area associated w/Lagrange multiplier vertex.
   */
  void _constrainSolnSpace1D(double_array* dLagrangeTpdt,
           const double_array& slip,
           const double_array& slipRate,
           const double_array& tractionTpdt,
           const double area);

  /** Constrain solution space with lumped Jacobian in 2-D.
   *
   * @param dLagrangeTpdt Adjustment to Lagrange multiplier.
   * @param slip Slip assoc. w/Lagrange multiplier vertex.
   * @param slipRate Slip rate assoc. w/Lagrange multiplier vertex.
   * @param tractionTpdt Fault traction assoc. w/Lagrange multiplier vertex.
   * @param area Fault area associated w/Lagrange multiplier vertex.
   */
  void _constrainSolnSpace2D(double_array* dLagrangeTpdt,
           const double_array& slip,
           const double_array& slipRate,
           const double_array& tractionTpdt,
           const double area);

  /** Constrain solution space with lumped Jacobian in 3-D.
   *
   * @param dLagrangeTpdt Adjustment to Lagrange multiplier.
   * @param slip Slip assoc. w/Lagrange multiplier vertex.
   * @param slipRate Slip rate assoc. w/Lagrange multiplier vertex.
   * @param tractionTpdt Fault traction assoc. w/Lagrange multiplier vertex.
   * @param area Fault area associated w/Lagrange multiplier vertex.
   */
  void _constrainSolnSpace3D(double_array* dLagrangeTpdt,
           const double_array& slip,
           const double_array& slipRate,
           const double_array& tractionTpdt,
           const double area);

  // PRIVATE MEMBERS ////////////////////////////////////////////////////
private :

  /// Database for initial tractions.
  spatialdata::spatialdb::SpatialDB* _dbInitialTract;

  /// To identify constitutive model
  friction::FrictionModel* _friction;

  /// Sparse matrix for sensitivity solve.
  topology::Jacobian* _jacobian;

  PetscKSP _ksp; ///< PETSc KSP linear solver for sensitivity problem.

// NOT IMPLEMENTED ////////////////////////////////////////////////////
private :

  /// Not implemented
  FaultCohesiveDyn(const FaultCohesiveDyn&);

  /// Not implemented
  const FaultCohesiveDyn& operator=(const FaultCohesiveDyn&);

}; // class FaultCohesiveDyn

#endif // pylith_faults_faultcohesivedyn_hh


// End of file 
