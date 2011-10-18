// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2011 University of California, Davis
//
// See COPYING for license information.
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
   *   be up-dip direction; applies to fault surfaces in 2-D and 3-D).
   */
  void initialize(const topology::Mesh& mesh,
		  const PylithScalar upDir[3]);

  /** Integrate contributions to residual term (r) for operator that
   * do not require assembly across processors.
   *
   * Initial tractions (if specified) are already assembled and
   * contribute to the residual like Neumann boundary conditions.
   *
   * @param residual Field containing values for residual
   * @param t Current time
   * @param fields Solution fields
   */
  virtual
  void integrateResidual(const topology::Field<topology::Mesh>& residual,
			 const PylithScalar t,
			 topology::SolutionFields* const fields);

  /** Update state variables as needed.
   *
   * @param t Current time
   * @param fields Solution fields
   * @param mesh Finite-element mesh
   */
  void updateStateVars(const PylithScalar t,
		       topology::SolutionFields* const fields);

  /** Constrain solution space based on friction.
   *
   * @param fields Solution fields.
   * @param t Current time.
   * @param jacobian Sparse matrix for system Jacobian.
   */
  void constrainSolnSpace(topology::SolutionFields* const fields,
			  const PylithScalar t,
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

  /** Update relative velocity associated with Lagrange vertex k
   * corresponding to diffential velocity between conventional
   * vertices i and j.
   *
   * @param fields Solution fields.
   */
  void _updateVelRel(const topology::SolutionFields& fields);

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
  void _sensitivitySolveLumped1D(scalar_array* slip,
                                 const scalar_array& dLagrangeTpdt,
                                 const scalar_array& jacobianN,
                                 const scalar_array& jacobianP);

  /** Solve slip/Lagrange multiplier sensitivity problem for case of lumped Jacobian in 2-D.
   *
   * @param slip Adjustment to slip assoc. w/Lagrange multiplier vertex.
   * @param dLagrangeTpdt Change in Lagrange multiplier.
   * @param jacobianN Jacobian for vertex on - side of the fault.
   * @param jacobianP Jacobian for vertex on + side of the fault.
   */
  void _sensitivitySolveLumped2D(scalar_array* slip,
                                 const scalar_array& dLagrangeTpdt,
                                 const scalar_array& jacobianN,
                                 const scalar_array& jacobianP);

  /** Solve slip/Lagrange multiplier sensitivity problem for case of lumped Jacobian in 3-D.
   *
   * @param slip Adjustment to slip assoc. w/Lagrange multiplier vertex.
   * @param dLagrangeTpdt Change in Lagrange multiplier.
   * @param jacobianN Jacobian for vertex on - side of the fault.
   * @param jacobianP Jacobian for vertex on + side of the fault.
   */
  void _sensitivitySolveLumped3D(scalar_array* slip,
                                 const scalar_array& dLagrangeTpdt,
                                 const scalar_array& jacobianN,
                                 const scalar_array& jacobianP);

  /** Constrain solution space with lumped Jacobian in 1-D.
   *
   * @param dLagrangeTpdt Adjustment to Lagrange multiplier.
   * @param slip Slip assoc. w/Lagrange multiplier vertex.
   * @param slipRate Slip rate assoc. w/Lagrange multiplier vertex.
   * @param tractionTpdt Fault traction assoc. w/Lagrange multiplier vertex.
   */
  void _constrainSolnSpace1D(scalar_array* dLagrangeTpdt,
           const scalar_array& slip,
           const scalar_array& slipRate,
           const scalar_array& tractionTpdt);

  /** Constrain solution space with lumped Jacobian in 2-D.
   *
   * @param dLagrangeTpdt Adjustment to Lagrange multiplier.
   * @param slip Slip assoc. w/Lagrange multiplier vertex.
   * @param slipRate Slip rate assoc. w/Lagrange multiplier vertex.
   * @param tractionTpdt Fault traction assoc. w/Lagrange multiplier vertex.
   */
  void _constrainSolnSpace2D(scalar_array* dLagrangeTpdt,
           const scalar_array& slip,
           const scalar_array& slipRate,
           const scalar_array& tractionTpdt);

  /** Constrain solution space with lumped Jacobian in 3-D.
   *
   * @param dLagrangeTpdt Adjustment to Lagrange multiplier.
   * @param slip Slip assoc. w/Lagrange multiplier vertex.
   * @param slipRate Slip rate assoc. w/Lagrange multiplier vertex.
   * @param tractionTpdt Fault traction assoc. w/Lagrange multiplier vertex.
   */
  void _constrainSolnSpace3D(scalar_array* dLagrangeTpdt,
           const scalar_array& slip,
           const scalar_array& slipRate,
           const scalar_array& tractionTpdt);

  // PRIVATE MEMBERS ////////////////////////////////////////////////////
private :

  /// Database for initial tractions.
  spatialdata::spatialdb::SpatialDB* _dbInitialTract;

  /// To identify constitutive model
  friction::FrictionModel* _friction;

  /// Sparse matrix for sensitivity solve.
  topology::Jacobian* _jacobian;

  PetscKSP _ksp; ///< PETSc KSP linear solver for sensitivity problem.

  /// Minimum resolvable value accounting for roundoff errors
  static const PylithScalar _zeroTolerance;

// NOT IMPLEMENTED ////////////////////////////////////////////////////
private :

  /// Not implemented
  FaultCohesiveDyn(const FaultCohesiveDyn&);

  /// Not implemented
  const FaultCohesiveDyn& operator=(const FaultCohesiveDyn&);

}; // class FaultCohesiveDyn

#endif // pylith_faults_faultcohesivedyn_hh


// End of file 
