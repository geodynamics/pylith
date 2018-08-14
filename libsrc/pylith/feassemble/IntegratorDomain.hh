// -*- C++ -*-
//
// ======================================================================
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, Rice University
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2018 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

/**
 * @file libsrc/feassemble/IntegratorDomain.hh
 *
 * @brief Object for finite-element integration over a subset (material) of the simulation domain.
 */

#if !defined(pylith_feassemble_integratordomain_hh)
#define pylith_feassemble_integratordomain_hh

#include "feassemblefwd.hh" // forward declarations

#include "pylith/feassemble/Integrator.hh" // ISA Integrator

class pylith::feassemble::IntegratorDomain : public pylith::feassemble::Integrator {
    friend class TestIntegratorDomain; // unit testing

    // PUBLIC MEMBERS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    /// Constructor
    IntegratorDomain(const pylith::problems::Physics* physics);

    /// Destructor
    virtual ~IntegratorDomain(void);

    /// Deallocate PETSc and local data structures.
    virtual
    void deallocate(void);

    /** Get spatial dimension of material.
     *
     * @returns Spatial dimension.
     */
    int getDimension(void) const;

    /** Set value of label material-id used to identify material cells.
     *
     * @param value Material identifier
     */
    void setMaterialId(const int value);

    /** Get value of label material-id used to identify material cells.
     *
     * @returns Material identifier
     */
    int getMaterialId(void) const;

    /** Set gravity field.
     *
     * @param g Gravity field.
     */
    void setGravityField(spatialdata::spatialdb::GravityField* const g);

    /** Add kernels for RHS residual.
     *
     * @param kernels Array of kernerls for computing the RHS residual.
     */
    void setKernelsRHSResidual(const std::vector<ResidualKernels>& kernels);

    /** Add kernels for RHS Jacobian.
     *
     * @param kernels Array of kernerls for computing the RHS Jacobian.
     */
    void setKernelsRHSJacobian(const std::vector<JacobianKernels>& kernels);

    /** Add kernels for RHS residual.
     *
     * @param kernels Array of kernerls for computing the RHS residual.
     */
    void setKernelsLHSResidual(const std::vector<ResidualKernels>& kernels);

    /** Add kernels for LHS Jacobian.
     *
     * @param kernels Array of kernerls for computing the LHS Jacobian.
     */
    void setKernelsLHSJacobian(const std::vector<JacobianKernels>& kernels);

    /** Add kernels for updating state variables.
     *
     * @param kernels Array of kernels for updating state variables.
     */
    void setKernelsUpdateStateVars(const std::vector<ProjectKernels>& kernels);

    /** Add kernels for computing derived field.
     *
     * @param kernels Array of kernels for computing derived field.
     */
    void setKernelsDerivedField(const std::vector<ProjectKernels>& kernels);

    /** Initialize integrator.
     *
     * @param[in] solution Solution field (layout).
     */
    void initialize(const pylith::topology::Field& solution);

    /** Compute RHS residual for G(t,s).
     *
     * @param[out] residual Field for residual.
     * @param[in] t Current time.
     * @param[in] dt Current time step.
     * @param[in] solution Field with current trial solution.
     */
    void computeRHSResidual(pylith::topology::Field* residual,
                            const PylithReal t,
                            const PylithReal dt,
                            const pylith::topology::Field& solution);

    /** Compute RHS Jacobian and preconditioner for G(t,s).
     *
     * @param[out] jacobianMat PETSc Mat with Jacobian sparse matrix.
     * @param[out] precondMat PETSc Mat with Jacobian preconditioning sparse matrix.
     * @param[in] t Current time.
     * @param[in] dt Current time step.
     * @param[in] solution Field with current trial solution.
     */
    void computeRHSJacobian(PetscMat jacobianMat,
                            PetscMat preconMat,
                            const PylithReal t,
                            const PylithReal dt,
                            const pylith::topology::Field& solution);

    /** Compute LHS residual for F(t,s,\dot{s}).
     *
     * @param[out] residual Field for residual.
     * @param[in] t Current time.
     * @param[in] dt Current time step.
     * @param[in] solution Field with current trial solution.
     * @param[in] solutionDot Field with time derivative of current trial solution.
     */
    void computeLHSResidual(pylith::topology::Field* residual,
                            const PylithReal t,
                            const PylithReal dt,
                            const pylith::topology::Field& solution,
                            const pylith::topology::Field& solutionDot);

    /** Compute LHS Jacobian and preconditioner for F(t,s,\dot{s}) with implicit time-stepping.
     *
     * @param[out] jacobianMat PETSc Mat with Jacobian sparse matrix.
     * @param[out] precondMat PETSc Mat with Jacobian preconditioning sparse matrix.
     * @param[in] t Current time.
     * @param[in] dt Current time step.
     * @param[in] s_tshift Scale for time derivative.
     * @param[in] solution Field with current trial solution.
     * @param[in] solutionDot Field with time derivative of current trial solution.
     */
    void computeLHSJacobian(PetscMat jacobianMat,
                            PetscMat precondMat,
                            const PylithReal t,
                            const PylithReal dt,
                            const PylithReal s_tshift,
                            const pylith::topology::Field& solution,
                            const pylith::topology::Field& solutionDot);

    /** Compute inverse of lumped LHS Jacobian for F(t,s,\dot{s}) with explicit time-stepping.
     *
     * @param[out] jacobianInv Inverse of lumped Jacobian as a field.
     * @param[in] t Current time.
     * @param[in] dt Current time step.
     * @param[in] s_tshift Scale for time derivative.
     * @param[in] solution Field with current trial solution.
     */
    void computeLHSJacobianLumpedInv(pylith::topology::Field* jacobianInv,
                                     const PylithReal t,
                                     const PylithReal dt,
                                     const PylithReal s_tshift,
                                     const pylith::topology::Field& solution);

    // PROTECTED METHODS ///////////////////////////////////////////////////////////////////////////////////////////////
protected:

    /** Compute residual using current kernels.
     *
     * @param[out] residual Field for residual.
     * @param[in] t Current time.
     * @param[in] dt Current time step.
     * @param[in] solution Field with current trial solution.
     * @param[in] solutionDot Field with time derivative of current trial solution.
     */
    void _computeResidual(pylith::topology::Field* residual,
                          const PylithReal t,
                          const PylithReal dt,
                          const pylith::topology::Field& solution,
                          const pylith::topology::Field& solutionDot);

    /** Compute Jacobian using current kernels.
     *
     * @param[out] jacobianMat PETSc Mat with Jacobian sparse matrix.
     * @param[out] precondMat PETSc Mat with Jacobian preconditioning sparse matrix.
     * @param[in] t Current time.
     * @param[in] dt Current time step.
     * @param[in] s_tshift Scale for time derivative.
     * @param[in] solution Field with current trial solution.
     * @param[in] solutionDot Field with time derivative of current trial solution.
     */
    void _computeJacobian(PetscMat jacobianMat,
                          PetscMat precondMat,
                          const PylithReal t,
                          const PylithReal dt,
                          const PylithReal s_tshift,
                          const pylith::topology::Field& solution,
                          const pylith::topology::Field& solutionDot);

    /** Update state variables as needed.
     *
     * @param[in] t Current time.
     * @param[in] dt Current time step.
     * @param[in] solution Field with current trial solution.
     */
    void _updateStateVars(const PylithReal t,
                          const PylithReal dt,
                          const pylith::topology::Field& solution);

    /** Compute fields derived from solution and auxiliary field.
     *
     * @param[in] t Current time.
     * @param[in] dt Current time step.
     * @param[in] solution Field with current trial solution.
     */
    void _computeDerivedFields(const PylithReal t,
                               const PylithReal dt,
                               const pylith::topology::Field& solution);

    // PRIVATE MEMBERS /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    std::vector<ResidualKernels> _kernelsRHSResidual; ///< kernels for RHS residual.
    std::vector<ResidualKernels> _kernelsRHSResidual; ///< kernels for LHS residual.

    std::vector<JacobianKernels> _kernelsRHSJacobian; ///< kernels for RHS Jacobian.
    std::vector<JacobianKernels> _kernelsLHSJacobian; ///> kernels for LHS Jacobian.

    std::vector<ProjectKernels> _kernelsUpstateStateVars; ///< kernels for updating state variables.
    std::vector<ProjectKernels> _kernelsDerivedField; ///< kernels for computing derived field.

    const int _dimension;
    int _materialId;

    // NOT IMPLEMENTED /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    IntegratorDomain(const IntegratorDomain&); ///< Not implemented.
    const IntegratorDomain& operator=(const IntegratorDomain&); ///< Not implemented.

}; // IntegratorDomain

#endif // pylith_feassemble_integratordomain_hh

// End of file
