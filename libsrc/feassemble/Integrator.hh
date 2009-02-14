// -*- C++ -*-
//
// ======================================================================
//
//                           Brad T. Aagaard
//                        U.S. Geological Survey
//
// {LicenseText}
//
// ======================================================================
//

/**
 * @file pylith/feassemble/Integrator.hh
 *
 * @brief Abstract base class for integration of finite-element
 * actions.
 *
 * Note: Each object operates on a single finite-element family, which
 * is defined by the quadrature and a database of material property
 * parameters.
 */

#if !defined(pylith_feassemble_integrator_hh)
#define pylith_feassemble_integrator_hh

// Include directives ---------------------------------------------------
#include "feassemblefwd.hh" // forward declarations

#include "pylith/topology/topologyfwd.hh" // USES Mesh, Field, SolutionFields
#include "spatialdata/spatialdb/spatialdbfwd.hh" // USES GravityField
#include "spatialdata/units/unitsfwd.hh" // USES Nondimensional

#include "pylith/utils/petscfwd.h" // USES PetscMat
#include "pylith/utils/array.hh" // HASA double_array

// Integrator -----------------------------------------------------------
template<typename quadrature_type>
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

  /** Set quadrature for integrating finite-element
   * quantities. Quadrature should already be initialized.
   *
   * @param q Quadrature for integrating.
   */
  void quadrature(const quadrature_type* q);

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
  void timeStep(const double dt);

  /** Get stable time step for advancing from time t to time t+dt.
   *
   * Default is MAXFLOAT (or 1.0e+30 if MAXFLOAT is not defined in math.h).
   *
   * @returns Time step
   */
  virtual
  double stableTimeStep(void) const;

  /** Check whether Jacobian needs to be recomputed.
   *
   * @returns True if Jacobian needs to be recomputed, false otherwise.
   */
  virtual
  bool needNewJacobian(void) const;

  /** Set flag for setting constraints for total field solution or
   *  incremental field solution.
   *
   * @param flag True if using incremental solution, false otherwise.
   */
  virtual
  void useSolnIncr(const bool flag);

  /** Integrate contributions to residual term (r) for operator.
   *
   * @param residual Field containing values for residual
   * @param t Current time
   * @param fields Solution fields
   */
  virtual 
  void integrateResidual(const topology::Field<topology::Mesh>& residual,
			 const double t,
			 topology::SolutionFields* const fields);

  /** Integrate contributions to Jacobian matrix (A) associated with
   * operator.
   *
   * @param jacobian Sparse matrix for Jacobian of system.
   * @param t Current time
   * @param fields Solution fields
   */
  virtual
  void integrateJacobian(PetscMat* jacobian,
			 const double t,
			 topology::SolutionFields* const fields);

  /** Integrate contributions to residual term (r) for operator that
   * do not require assembly over cells, vertices, or processors.
   *
   * @param residual Field containing values for residual
   * @param t Current time
   * @param fields Solution fields
   */
  virtual 
  void integrateResidualAssembled(const topology::Field<topology::Mesh>& residual,
				  const double t,
				  topology::SolutionFields* const fields);

  /** Integrate contributions to Jacobian matrix (A) associated with
   * operator that do not require assembly over cells, vertices, or
   * processors
   *
   * @param jacobian Sparse matrix for Jacobian of system.
   * @param t Current time
   * @param fields Solution fields
   */
  virtual
  void integrateJacobianAssembled(PetscMat* jacobian,
				  const double t,
				  topology::SolutionFields* const fields);

  /** Update state variables as needed.
   *
   * @param t Current time
   * @param fields Solution fields
   * @param mesh Finite-element mesh
   */
  virtual
  void updateState(const double t,
		   topology::SolutionFields* const fields);

  /** Verify configuration is acceptable.
   *
   * @param mesh Finite-element mesh
   */
  virtual
  void verifyConfiguration(const topology::Mesh& mesh) const = 0;

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

// PROTECTED MEMBERS ////////////////////////////////////////////////////
protected :

  double _dt; ///< Time step for t -> t+dt

  quadrature_type* _quadrature; ///< Quadrature for integrating finite-element

  spatialdata::units::Nondimensional* _normalizer; ///< Nondimensionalizer.
  spatialdata::spatialdb::GravityField* _gravityField; ///< Gravity field.

  /// Vector local to cell containing result of integration action
  double_array _cellVector;

  /// Matrix local to cell containing result of integration
  double_array _cellMatrix;

  /// True if we need to recompute Jacobian for operator, false otherwise.
  /// Default is false;
  bool _needNewJacobian;

  /// Flag indicating whether to set constraints for a total field
  /// solution or an incremental field solution
  bool _useSolnIncr;

// NOT IMPLEMENTED //////////////////////////////////////////////////////
private :

  Integrator(const Integrator& i); ///< Not implemented
  const Integrator& operator=(const Integrator&); ///< Not implemented

}; // Integrator

#include "Integrator.icc" // inline methods
#include "Integrator.cc" // template methods

#endif // pylith_feassemble_integrator_hh


// End of file 
