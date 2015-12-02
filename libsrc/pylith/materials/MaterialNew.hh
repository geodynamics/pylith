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
// Copyright (c) 2010-2015 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

/** @file libsrc/materials/MaterialNew.hh
 *
 * @brief C++ abstract base class for Material object.
 */

#if !defined(pylith_materials_materialnew_hh)
#define pylith_materials_materialnew_hh

// Include directives ---------------------------------------------------
#include "materialsfwd.hh" // forward declarations

#include "pylith/feassemble/IntegratorPointwise.hh" // ISA IntegratorPointwise

#include "pylith/topology/FieldBase.hh" // HASA FieldBase::DiscretizeInfo

#include "spatialdata/spatialdb/spatialdbfwd.hh" // forward declarations

#include <string> // HASA std::string

// Material -------------------------------------------------------------
/** @brief C++ abstract base class for Material object.
 *
 * Interface definition for a material. A material encapsulates both
 * the rheology as well as the governing equation.
 *
 * An individual material must abide by specific rules for the
 * interface, especially the order of the fields in the solution.
 *
 * Elasticity:
 *   + displacement, [velocity, Lagrange multipliers]
 *
 * Incompressible elasticity
 *   + displacement, pressure, [velocity, Lagrange multipliers]
 */

class pylith::materials::MaterialNew : public pylith::feassemble::IntegratorPointwise
{ // class Material
  friend class TestMaterialNew; // unit testing

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /** Default constructor.
   *
   * @param dimension Spatial dimension associated with material.
   */
  MaterialNew(const int dimension);

  /// Destructor.
  virtual
  ~MaterialNew(void);

  /// Deallocate PETSc and local data structures.
  virtual
  void deallocate(void);
  
  /** Get spatial dimension of material.
   *
   * @returns Spatial dimension.
   */
  int dimension(void) const;

  /** Set identifier of material.
   *
   * @param value Material identifier
   */
  void id(const int value);

  /** Get identifier of material.
   *
   * @returns Material identifier
   */
  int id(void) const;

  /** Set label of material.
   *
   * @param value Label of material
   */
  void label(const char* value);

  /** Get label of material.
   *
   * @returns Label of material
   */
  const char* label(void) const;

  /** Initialize material. Setup auxiliary fields.
   *
   * @param solution Solution field.
   */
  void initialize(const pylith::topology::Field& solution);
  
  /** Compute RHS residual for G(t,s).
   *
   * @param[out] residual Residual field.
   * @param[in] t Current time.
   * @param[in] dt Current time step.
   * @param[in] solution Current trial solution.
   */
  void computeRHSResidual(pylith::topology::Field* residual,
			  const PylithReal t,
			  const PylithReal dt,
			  const pylith::topology::Field& solution);

  /** Compute RHS Jacobian and preconditioner for G(t,s).
   *
   * @param[out] jacobian Jacobian sparse matrix.
   * @param[out] preconditioner Jacobian preconditioning sparse matrix.
   * @param[in] t Current time.
   * @param[in] dt Current time step.
   * @param[in] solution Current trial solution.
   */
  void computeRHSJacobian(pylith::topology::Jacobian* jacobian,
			  pylith::topology::Jacobian* preconditioner,
			  const PylithReal t,
			  const PylithReal dt,
			  const pylith::topology::Field& solution);

  /** Compute LHS residual for F(t,s,\dot{s}).
   *
   * @param[out] residual Residual field.
   * @param[in] t Current time.
   * @param[in] dt Current time step.
   * @param[in] solution Current trial solution.
   */
  void computeLHSResidual(pylith::topology::Field* residual,
			  const PylithReal t,
			  const PylithReal dt,
			  const pylith::topology::Field& solution);

  /** Compute LHS Jacobian for F(t,s,\dot{s}) with implicit time-stepping.
   *
   * @param[out] jacobian Jacobian sparse matrix.
   * @param[out] preconditioner Jacobian preconditioning sparse matrix.
   * @param[in] t Current time.
   * @param[in] dt Current time step.
   * @param[in] solution Current trial solution.
   */
  void computeLHSJacobianImplicit(pylith::topology::Jacobian* jacobian,
				  pylith::topology::Jacobian* preconditioner,
				  const PylithReal t,
				  const PylithReal dt,
				  const pylith::topology::Field& solution);


  /** Compute LHS Jacobian for F(t,s,\dot{s}) with explicit time-stepping.
   *
   * @param[out] jacobian Jacobian sparse matrix.
   * @param[out] preconditioner Jacobian preconditioning sparse matrix.
   * @param[in] t Current time.
   * @param[in] dt Current time step.
   * @param[in] solution Current trial solution.
   */
  void computeLHSJacobianExplicit(pylith::topology::Jacobian* jacobian,
				  pylith::topology::Jacobian* preconditioner,
				  const PylithReal t,
				  const PylithReal dt,
				  const pylith::topology::Field& solution);


  /** Update state variables as needed.
   *
   * @param[in] solution Solution field.
   */
  void updateStateVars(const pylith::topology::Field& solution);

  // PROTECTED METHODS //////////////////////////////////////////////////
protected :

  /* Compute residual using current kernels.
   *
   * @param[out] residual Residual field.
   * @param[in] t Current time.
   * @param[in] dt Current time step.
   * @param[in] solution Current trial solution.
   */
  void _computeResidual(pylith::topology::Field* residual,
			const PylithReal t,
			const PylithReal dt,
			const pylith::topology::Field& solution);
  
  /* Compute Jacobian using current kernels.
   *
   * @param[out] jacobian Jacobian sparse matrix.
   * @param[out] preconditioner Jacobian preconditioning sparse matrix.
   * @param[in] t Current time.
   * @param[in] dt Current time step.
   * @param[in] solution Current trial solution.
   */
  void _computeJacobian(pylith::topology::Jacobian* jacobian,
			pylith::topology::Jacobian* preconditioner,
			const PylithReal t,
			const PylithReal dt,
			const pylith::topology::Field& solution);
  
  /** Setup auxiliary subfields (discretization and query fns).
   * 
   * Create subfields in auxiliary fields (includes name of the field,
   * vector field type, discretization, and scale for
   * nondimensionalization) and set query functions for filling them
   * from a spatial database.
   *
   * @attention The order of the calls to subfieldAdd() must match the
   * order of the auxiliary fields in the FE kernels.
   */
  virtual
  void _auxFieldsSetup(void) = 0;

  /** Set kernels for RHS residual G(t,s).
   *
   * Potentially, there are g0 and g1 kernels for each equation. If no
   * kernel is needed, then set the kernel function to NULL.
   *
   * @param solution Solution field.
   */
  virtual
  void _setFEKernelsRHSResidual(const topology::Field& solution) const = 0;


  /** Set kernels for RHS Jacobian G(t,s).
   *
   * Potentially, there are Jg0, Jg1, Jg2, and Jg3 kernels for each
   * combination of equations. If no kernel is needed, then set the
   * kernel function to NULL.
   *
   * - Jg0(ifield, jfield)
   * - Jg1(ifield, jfield, jdim)
   * - Jg2(ifield, jfield, idim)
   * - Jg3(ifield, jfield, idim, jdim)
   *
   * @param solution Solution field.
   */
  virtual
  void _setFEKernelsRHSJacobian(const topology::Field& solution) const = 0;


  /** Set kernels for LHS residual F(t,s,\dot{s}).
   *
   * Potentially, there are f0 and f1 kernels for each equation. If no
   * kernel is needed, then set the kernel function to NULL.
   *
   * @param solution Solution field.
   */
  virtual
  void _setFEKernelsLHSResidual(const topology::Field& solution) const = 0;


  /** Set kernels for LHS Jacobian F(t,s,\dot{s}) when implicit time-stepping.
   *
   * - Jf0(ifield, jfield)
   * - Jf1(ifield, jfield, jdim)
   * - Jf2(ifield, jfield, idim)
   * - Jf3(ifield, jfield, idim, jdim)
   *
   * @param solution Solution field.
   */
  virtual
  void _setFEKernelsLHSJacobianImplicit(const topology::Field& solution) const = 0;


  /** Set kernels for LHS Jacobian F(t,s,\dot{s}) when explicit time-stepping.
   *
   * - Jf0(ifield, jfield)
   * - Jf1(ifield, jfield, jdim)
   * - Jf2(ifield, jfield, idim)
   * - Jf3(ifield, jfield, idim, jdim)
   *
   * @param solution Solution field.
   */
  virtual
  void _setFEKernelsLHSJacobianExplicit(const topology::Field& solution) const = 0;


  // PROTECTED MEMBERS //////////////////////////////////////////////////
protected :

  topology::StratumIS* _materialIS; ///< Index set for material cells.

  // PRIVATE MEMBERS ////////////////////////////////////////////////////
private :

  const int _dimension; ///< Spatial dimension of material.
  int _id; ///< Material identifier.
  std::string _label; ///< Label of material.

  // NOT IMPLEMENTED ////////////////////////////////////////////////////
private :

  MaterialNew(const MaterialNew&); ///< Not implemented.
  const MaterialNew& operator=(const MaterialNew&); ///< Not implemented

}; // class Material

#include "MaterialNew.icc" // inline methods

#endif // pylith_materials_material_hh


// End of file 
