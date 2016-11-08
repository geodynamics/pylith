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
// Copyright (c) 2010-2016 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

/** @file modulesrc/feassemble/IntegratorPointwise.i
 *
 * @brief Python interface to C++ abstract IntegratorPointwise object.
 */

namespace pylith {
namespace feassemble {

class IntegratorPointwise
{     // IntegratorPointwise

// PUBLIC MEMBERS ///////////////////////////////////////////////////////
public:

/// Constructor
IntegratorPointwise(void);

/// Destructor
virtual
~IntegratorPointwise(void);

/// Deallocate PETSc and local data structures.
virtual
void deallocate(void);

/** Get auxiliary fields.
 *
 * @return field Field over material.
 */
const pylith::topology::Field& auxFields(void) const;

/** Check whether material has a given auxiliary field.
 *
 * @param[in] name Name of field.
 *
 * @returns True if material has auxiliary field, false otherwise.
 */
bool hasAuxField(const char* name);

/** Get auxiliary field.
 *
 * @param[in] field Field over material.
 * @param[in] name Name of field to retrieve.
 */
  void getAuxField(pylith::topology::Field *field,
                 const char* name) const;

/** Set spatial database for auxiliary fields.
 *
 * @param[in] value Pointer to database.
 */
void auxFieldsDB(spatialdata::spatialdb::SpatialDB* value);

/** Set discretization information for auxiliary subfield.
 *
 * @param[in] name Name of auxiliary subfield.
 * @feInfo Discretization information for subfield.
 */
void auxFieldDiscretization(const char* name,
                            const pylith::topology::FieldBase::DiscretizeInfo& feInfo);

/** Get discretization information for auxiliary subfield.
 *
 * @param[in] name Name of subfield.
 * @return Discretization information for auxiliary subfield. If
 * discretization information was not set, then use "default".
 */
const pylith::topology::FieldBase::DiscretizeInfo& auxFieldDiscretization(const char* name) const;

/** Check whether RHS Jacobian needs to be recomputed.
 *
 * @returns True if Jacobian needs to be recomputed, false otherwise.
 */
virtual
bool needNewRHSJacobian(void) const;

/** Check whether LHS Jacobian needs to be recomputed.
 *
 * @returns True if Jacobian needs to be recomputed, false otherwise.
 */
virtual
bool needNewLHSJacobian(void) const;

/** Set manager of scales used to nondimensionalize problem.
 *
 * @param dim Nondimensionalizer.
 */
void normalizer(const spatialdata::units::Nondimensional& dim);

/** Verify configuration is acceptable.
 *
 * @param[in] mesh Finite-element mesh
 */
virtual
void verifyConfiguration(const pylith::topology::Mesh& mesh) const;

/** Verify constraints are acceptable.
 *
 * @param field Solution field.
 */
virtual
void checkConstraints(const pylith::topology::Field& solution) const;

/** Initialize integrator.
 *
 * @param[in] solution Solution field (layout).
 */
virtual
void initialize(const pylith::topology::Field& solution) = 0;

/** Compute RHS residual for G(t,s).
 *
 * @param[out] residualVec PETSc Vec for residual field.
 * @param[in] t Current time.
 * @param[in] dt Current time step.
 * @param[in] solution Field with current trial solution.
 */
virtual
void computeRHSResidual(PetscVec residualVec,
                        const PylithReal t,
                        const PylithReal dt,
                        const pylith::topology::Field& solution) = 0;

/** Compute RHS Jacobian and preconditioner for G(t,s).
 *
 * @param[out] jacobianMat PETSc Mat with Jacobian sparse matrix.
 * @param[out] precondMat PETSc Mat with Jacobian preconditioning sparse matrix.
 * @param[in] t Current time.
 * @param[in] dt Current time step.
 * @param[in] solution Field with current trial solution.
 */
virtual
void computeRHSJacobian(PetscMat jacobianMat,
                        PetscMat preconMat,
                        const PylithReal t,
                        const PylithReal dt,
                        const pylith::topology::Field& solution) = 0;

/** Compute LHS residual for F(t,s,\dot{s}).
 *
 * @param[out] residualVec PETSc Vec for residual field.
 * @param[in] t Current time.
 * @param[in] dt Current time step.
 * @param[in] solution Field with current trial solution.
 * @param[in] solutionDotVec PETSc Vec with time derivative of current trial solution.
 */
virtual
void computeLHSResidual(PetscVec residualVec,
                        const PylithReal t,
                        const PylithReal dt,
                        const pylith::topology::Field& solution,
                        PetscVec solutionDotVec) = 0;

/** Compute LHS Jacobian and preconditioner for F(t,s,\dot{s}) with implicit time-stepping.
 *
 * @param[out] jacobianMat PETSc Mat with Jacobian sparse matrix.
 * @param[out] precondMat PETSc Mat with Jacobian preconditioning sparse matrix.
 * @param[in] t Current time.
 * @param[in] dt Current time step.
 * @param[in] tshift Scale for time derivative.
 * @param[in] solution Field with current trial solution.
 * @param[in] solutionDotVec PETSc Vec with time derivative of current trial solution.
 */
virtual
void computeLHSJacobianImplicit(PetscMat jacobianMat,
                                PetscMat precondMat,
                                const PylithReal t,
                                const PylithReal dt,
                                const PylithReal tshift,
                                const pylith::topology::Field& solution,
                                PetscVec solutionDotVec) = 0;


/** Compute inverse of lumped LHS Jacobian for F(t,s,\dot{s}) with explicit time-stepping.
 *
 * @param[out] jacobianInv Inverse of lumped Jacobian as a field.
 * @param[in] t Current time.
 * @param[in] dt Current time step.
 * @param[in] solution Field with current trial solution.
 */
virtual
void computeLHSJacobianLumpedInv(pylith::topology::Field* jacobianInv,
                                 const PylithReal t,
                                 const PylithReal dt,
                                 const pylith::topology::Field& solution) = 0;


/** Update state variables as needed.
 *
 * @param[in] solution Field with current trial solution.
 */
virtual
void updateStateVars(const pylith::topology::Field& solution);


};     // IntegratorPointwise

}   // feassemble
} // pylith


// End of file
