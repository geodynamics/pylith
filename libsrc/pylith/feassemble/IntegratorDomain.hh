// -*- C++ -*-
//
// ======================================================================
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University at Buffalo
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2021 University of California, Davis
//
// See LICENSE.md for license information.
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
#include "pylith/utils/arrayfwd.hh" // HASA std::vector

class pylith::feassemble::IntegratorDomain : public pylith::feassemble::Integrator {
    friend class TestIntegratorDomain; // unit testing

    // PUBLIC STRUCTS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    /// Kernels (point-wise functions) for residual.
    struct ResidualKernels {
        std::string subfield; ///< Name of subfield.
        ResidualPart part; ///< Residual part (LHS or RHS).
        PetscPointFunc r0; ///< f0 (RHS) or g0 (LHS) function.
        PetscPointFunc r1; ///< f1 (RHS) or g1 (LHS) function.

        ResidualKernels(void) :
            subfield(""),
            part(RESIDUAL_LHS),
            r0(NULL),
            r1(NULL) {}


        ResidualKernels(const char* subfieldValue,
                        const ResidualPart partValue,
                        PetscPointFunc r0Value,
                        PetscPointFunc r1Value) :
            subfield(subfieldValue),
            part(partValue),
            r0(r0Value),
            r1(r1Value) {}


    }; // ResidualKernels

    /// Kernels (point-wise functions) for Jacobian;
    struct JacobianKernels {
        std::string subfieldTrial; ///< Name of subfield associated with trial function (row in Jacobian).
        std::string subfieldBasis; ///< Name of subfield associated with basis function (column in Jacobian).
        JacobianPart part; ///< Jacobian part (LHS or LHS lumped inverse).
        PetscPointJac j0; ///< J0 function.
        PetscPointJac j1; ///< J1 function.
        PetscPointJac j2; ///< J2 function.
        PetscPointJac j3; ///< J3 function.

        JacobianKernels(void) :
            subfieldTrial(""),
            subfieldBasis(""),
            part(JACOBIAN_LHS),
            j0(NULL),
            j1(NULL),
            j2(NULL),
            j3(NULL) {}


        JacobianKernels(const char* subfieldTrialValue,
                        const char* subfieldBasisValue,
                        JacobianPart partValue,
                        PetscPointJac j0Value,
                        PetscPointJac j1Value,
                        PetscPointJac j2Value,
                        PetscPointJac j3Value) :
            subfieldTrial(subfieldTrialValue),
            subfieldBasis(subfieldBasisValue),
            part(partValue),
            j0(j0Value),
            j1(j1Value),
            j2(j2Value),
            j3(j3Value) {}


    }; // JacobianKernels

    /// Project kernels (pointwise functions) for updating state variables or computing derived fields.
    struct ProjectKernels {
        std::string subfield; ///< Name of subfield for function.
        PetscPointFunc f; ///< Point-wise function.

        ProjectKernels(void) :
            subfield(""),
            f(NULL) {}


        ProjectKernels(const char* subfieldValue,
                       PetscPointFunc fValue) :
            subfield(subfieldValue),
            f(fValue) {}


    }; // ProjectKernels

    // PUBLIC MEMBERS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    /// Constructor
    IntegratorDomain(pylith::problems::Physics* const physics);

    /// Destructor
    virtual ~IntegratorDomain(void);

    /// Deallocate PETSc and local data structures.
    virtual
    void deallocate(void);

    /** Get mesh associated with integrator domain.
     *
     * @returns Mesh associated with integrator domain.
     */
    const pylith::topology::Mesh& getPhysicsDomainMesh(void) const;

    /** Set kernels for residual.
     *
     * @param[in] kernels Array of kernerls for computing the residual.
     * @param[in] solution Solution field.
     */
    void setKernelsResidual(const std::vector<ResidualKernels>& kernels,
                            const pylith::topology::Field& solution);

    /** Set kernels for Jacobian.
     *
     * @param[in] kernels Array of kernerls for computing the Jacobian.
     * @param[in] solution Solution field.
     */
    void setKernelsJacobian(const std::vector<JacobianKernels>& kernels,
                            const pylith::topology::Field& solution);

    /** Set kernels for updating state variables.
     *
     * @param kernels Array of kernels for updating state variables.
     */
    void setKernelsUpdateStateVars(const std::vector<ProjectKernels>& kernels);

    /** Set kernels for computing derived field.
     *
     * @param kernels Array of kernels for computing derived field.
     */
    void setKernelsDerivedField(const std::vector<ProjectKernels>& kernels);

    /** Initialize integration domain, auxiliary field, and derived field. Update observers.
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

    /** Update state variables as needed.
     *
     * @param[in] t Current time.
     * @param[in] dt Current time step.
     * @param[in] solution Field with current trial solution.
     */
    void _updateStateVars(const PylithReal t,
                          const PylithReal dt,
                          const pylith::topology::Field& solution);

    /** Compute field derived from solution and auxiliary field.
     *
     * @param[in] t Current time.
     * @param[in] dt Current time step.
     * @param[in] solution Field with current trial solution.
     */
    void _computeDerivedField(const PylithReal t,
                              const PylithReal dt,
                              const pylith::topology::Field& solution);

    // PRIVATE MEMBERS /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    std::vector<ProjectKernels> _kernelsUpdateStateVars; ///< kernels for updating state variables.
    std::vector<ProjectKernels> _kernelsDerivedField; ///< kernels for computing derived field.

    pylith::topology::Mesh* _materialMesh; ///< Mesh associated with material.

    pylith::feassemble::UpdateStateVars* _updateState; ///< Data structure for layout needed to update state vars.

    // NOT IMPLEMENTED /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    IntegratorDomain(void); ///< Not implemented.
    IntegratorDomain(const IntegratorDomain&); ///< Not implemented.
    const IntegratorDomain& operator=(const IntegratorDomain&); ///< Not implemented.

}; // IntegratorDomain

#endif // pylith_feassemble_integratordomain_hh

// End of file
