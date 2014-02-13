// -*- C++ -*-
//
// ======================================================================
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2014 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

/**
 * @file libsrc/feassemble/Integrator.hh
 *
 * @brief Abstract base class for integration of finite-element
 * actions.
 */

#if !defined(pylith_feassemble_integrator_hh)
#define pylith_feassemble_integrator_hh

// Include directives ---------------------------------------------------
#include "feassemblefwd.hh" // forward declarations

#include "pylith/topology/topologyfwd.hh" // USES Mesh, Field, SolutionFields
#include "pylith/utils/utilsfwd.hh" // HOLDSA EventLogger
#include "pylith/utils/petscfwd.h" // USES PetscMat

#include "spatialdata/spatialdb/spatialdbfwd.hh" // USES GravityField
#include "spatialdata/units/unitsfwd.hh" // USES Nondimensional

#include "pylith/utils/array.hh" // HASA scalar_array

// Integrator -----------------------------------------------------------
/** @brief Abstract base class for integration of finite-element
 * actions.
 *
 * Note: Each object operates on a single finite-element family, which
 * is defined by the quadrature and a database of material property
 * parameters.
 */
class pylith::feassemble::Integrator
{ // Integrator
  friend class TestIntegrator; // unit testing

// PUBLIC MEMBERS ///////////////////////////////////////////////////////
public :

  /// Constructor
  Integrator(void);

  /// Destructor
  virtual
  ~Integrator(void);

  /// Deallocate PETSc and local data structures.
  virtual
  void deallocate(void);
  
  /** Get the quadrature for integrating finite-element
   * quantities.
   *
   * @returns Quadrature for integrating.
   */
  const Quadrature& quadrature();

  /** Set quadrature for integrating finite-element
   * quantities. Quadrature should already be initialized.
   *
   * @param q Quadrature for integrating.
   */
  void quadrature(const Quadrature* q);

  /** Set manager of scales used to nondimensionalize problem.
   *
   * @param dim Nondimensionalizer.
   */
  void normalizer(const spatialdata::units::Nondimensional& dim);

  /** Set gravity field.
   *
   * @param g Gravity field.
   */
  void gravityField(spatialdata::spatialdb::GravityField* const gravityField);

  /** Set time step for advancing from time t to time t+dt.
   *
   * @param dt Time step
   */
  virtual
  void timeStep(const PylithScalar dt);

  /** Get stable time step for advancing from time t to time t+dt.
   *
   * Default is MAXFLOAT (or 1.0e+30 if MAXFLOAT is not defined in math.h).
   *
   * @param mesh Finite-element mesh.
   * @returns Time step
   */
  virtual
  PylithScalar stableTimeStep(const topology::Mesh& mesh);

  /** Check whether Jacobian needs to be recomputed.
   *
   * @returns True if Jacobian needs to be recomputed, false otherwise.
   */
  virtual
  bool needNewJacobian(void) const;

  /** Check whether integrator generates a symmetric Jacobian.
   *
   * @returns True if integrator generates symmetric Jacobian.
   */
  virtual
  bool isJacobianSymmetric(void) const;

  /** Initialize integrator.
   *
   * @param mesh Finite-element mesh.
   */
  virtual
  void initialize(const topology::Mesh& mesh);

  /** Setup DOF on solution field.
   *
   * @param field Solution field.
   */
  virtual
  void setupSolnDof(topology::Field* field);

  /** Integrate contributions to residual term (r) for operator.
   *
   * @param residual Field containing values for residual
   * @param t Current time
   * @param fields Solution fields
   */
  virtual 
  void integrateResidual(const topology::Field& residual,
			 const PylithScalar t,
			 topology::SolutionFields* const fields);

  /** Integrate contributions to residual term (r) for operator.
   *
   * @param residual Field containing values for residual
   * @param t Current time
   * @param fields Solution fields
   */
  virtual
  void integrateResidualLumped(const topology::Field& residual,
       const PylithScalar t,
       topology::SolutionFields* const fields);

  /** Integrate contributions to Jacobian matrix (A) associated with
   * operator.
   *
   * @param jacobian Sparse matrix for Jacobian of system.
   * @param t Current time
   * @param fields Solution fields
   */
  virtual
  void integrateJacobian(topology::Jacobian* jacobian,
			 const PylithScalar t,
			 topology::SolutionFields* const fields);

  /** Integrate contributions to Jacobian matrix (A) associated with
   * operator.
   *
   * @param jacobian Diagonal matrix (as field) for Jacobian of system.
   * @param t Current time
   * @param fields Solution fields
   */
  virtual
  void integrateJacobian(topology::Field* jacobian,
			 const PylithScalar t,
			 topology::SolutionFields* const fields);

  /** Integrate contributions to Jacobian matrix (A) associated with
   * operator.
   *
   * @param precondMatrix Custom preconditioning matrix.
   * @param jacobian Sparse matrix for Jacobian of system.
   * @param fields Solution fields
   */
  virtual
  void calcPreconditioner(PetscMat* const precondMatrix,
			  topology::Jacobian* const jacobian,
			  topology::SolutionFields* const fields);

  /** Update state variables as needed.
   *
   * @param t Current time
   * @param fields Solution fields
   * @param mesh Finite-element mesh
   */
  virtual
  void updateStateVars(const PylithScalar t,
		       topology::SolutionFields* const fields);

  /** Constrain solution space.
   *
   * @param fields Solution fields.
   * @param t Current time.
   * @param jacobian Sparse matrix for system Jacobian.
   */
  virtual
  void constrainSolnSpace(topology::SolutionFields* const fields,
			  const PylithScalar t,
			  const topology::Jacobian& jacobian);

  /** Adjust solution from solver with lumped Jacobian to match Lagrange
   *  multiplier constraints.
   *
   * @param fields Solution fields.
   * @param t Current time.
   * @param jacobian Jacobian of the system.
   */
  virtual
  void adjustSolnLumped(topology::SolutionFields* fields,
			const PylithScalar t,
			const topology::Field& jacobian);

  /** Verify configuration is acceptable.
   *
   * @param mesh Finite-element mesh
   */
  virtual
  void verifyConfiguration(const topology::Mesh& mesh) const = 0;

  /** Verify constraints are acceptable.
   *
   * @param field Solution field.
   */
  virtual
  void checkConstraints(const topology::Field& solution) const;

// PROTECTED METHODS ////////////////////////////////////////////////////
protected :

  /// Initialize vector containing result of integration action for cell.
  void _initCellVector(void);

  /// Zero out vector containing result of integration actions for cell.
  void _resetCellVector(void);

  /// Initialize matrix containing result of integration for cell.
  void _initCellMatrix(void);

  /// Zero out matrix containing result of integration for cell.
  void _resetCellMatrix(void);

  /// Lump cell matrix, putting the result in the cell vector using
  /// equivalent forces for rigid body motion.
  void _lumpCellMatrix(void);

// PROTECTED MEMBERS ////////////////////////////////////////////////////
protected :

  PylithScalar _dt; ///< Time step for t -> t+dt

  Quadrature* _quadrature; ///< Quadrature for integrating finite-element

  spatialdata::units::Nondimensional* _normalizer; ///< Nondimensionalizer.
  spatialdata::spatialdb::GravityField* _gravityField; ///< Gravity field.

  utils::EventLogger* _logger; ///< Event logger.

  /// Vector local to cell containing result of integration action
  scalar_array _cellVector;

  /// Matrix local to cell containing result of integration
  scalar_array _cellMatrix;

  /// True if we need to recompute Jacobian for operator, false otherwise.
  /// Default is false;
  bool _needNewJacobian;

  /// True if we need to compute velocity field, false otherwise.
  /// Default is false;
  bool _isJacobianSymmetric;

// NOT IMPLEMENTED //////////////////////////////////////////////////////
private :

  Integrator(const Integrator& i); ///< Not implemented
  const Integrator& operator=(const Integrator&); ///< Not implemented

}; // Integrator

#include "Integrator.icc" // inline methods

#endif // pylith_feassemble_integrator_hh


// End of file 
