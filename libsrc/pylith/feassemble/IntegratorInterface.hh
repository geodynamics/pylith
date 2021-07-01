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
 * @file libsrc/feassemble/IntegratorInterface.hh
 *
 * @brief Object for finite-element integration over an interior interface of the simulation domain.
 */

#if !defined(pylith_feassemble_integratorinterface_hh)
#define pylith_feassemble_integratorinterface_hh

#include "feassemblefwd.hh" // forward declarations

#include "pylith/feassemble/Integrator.hh" // ISA Integrator
#include "pylith/utils/arrayfwd.hh" // HASA std::vector

class pylith::feassemble::IntegratorInterface : public pylith::feassemble::Integrator {
    friend class TestIntegratorInterface; // unit testing

    // PUBLIC STRUCTS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    /// Kernels (point-wise functions) for residual.
    struct ResidualKernels {
        std::string subfield; ///< Name of subfield
        PetscBdPointFunc r0; ///< f0 (RHS) or g0 (LHS) function.
        PetscBdPointFunc r1; ///< f1 (RHS) or g1 (LHS) function.

        ResidualKernels(void) :
            subfield(""),
            r0(NULL),
            r1(NULL) {}


        ResidualKernels(const char* subfieldValue,
                        PetscBdPointFunc r0Value,
                        PetscBdPointFunc r1Value) :
            subfield(subfieldValue),
            r0(r0Value),
            r1(r1Value) {}


    }; // ResidualKernels

    /// Kernels (point-wise functions) for Jacobian;
    struct JacobianKernels {
        std::string subfieldTrial; ///< Name of subfield associated with trial function (row in Jacobian).
        std::string subfieldBasis; ///< Name of subfield associated with basis function (column in Jacobian).
        PetscBdPointJac j0; ///< J0 function.
        PetscBdPointJac j1; ///< J1 function.
        PetscBdPointJac j2; ///< J2 function.
        PetscBdPointJac j3; ///< J3 function.

        JacobianKernels(void) :
            subfieldTrial(""),
            subfieldBasis(""),
            j0(NULL),
            j1(NULL),
            j2(NULL),
            j3(NULL) {}


        JacobianKernels(const char* subfieldTrialValue,
                        const char* subfieldBasisValue,
                        PetscBdPointJac j0Value,
                        PetscBdPointJac j1Value,
                        PetscBdPointJac j2Value,
                        PetscBdPointJac j3Value) :
            subfieldTrial(subfieldTrialValue),
            subfieldBasis(subfieldBasisValue),
            j0(j0Value),
            j1(j1Value),
            j2(j2Value),
            j3(j3Value) {}


    }; // JacobianKernels

    // PUBLIC MEMBERS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    /// Constructor
    IntegratorInterface(pylith::problems::Physics* const physics);

    /// Destructor
    virtual ~IntegratorInterface(void);

    /// Deallocate PETSc and local data structures.
    virtual
    void deallocate(void);

    /** Set label marking interior interface surface.
     *
     * @param[in] value Label of surface (from mesh generator).
     */
    void setSurfaceMarkerLabel(const char* value);

    /** Get label marking interior interface surface.
     *
     * @returns Label of surface (from mesh generator).
     */
    const char* getSurfaceMarkerLabel(void) const;

    /** Get mesh associated with integrator domain.
     *
     * @returns Mesh associated with integrator domain.
     */
    const pylith::topology::Mesh& getPhysicsDomainMesh(void) const;

    /** Set kernels for RHS residual.
     *
     * @param kernels Array of kernels for computing the RHS residual.
     */
    void setKernelsRHSResidual(const std::vector<ResidualKernels>& kernels);

    /** Set kernels for LHS residual.
     *
     * @param kernels Array of kernels for computing the LHS residual.
     */
    void setKernelsLHSResidual(const std::vector<ResidualKernels>& kernels);

    /** Set kernels for LHS Jacobian.
     *
     * @param kernels Array of kernels for computing the LHS Jacobian.
     */
    void setKernelsLHSJacobian(const std::vector<JacobianKernels>& kernels);

    /** Initialize integration domain, auxiliary field, and derived field. Update observers.
     *
     * @param[in] solution Solution field (layout).
     */
    void initialize(const pylith::topology::Field& solution);

    /** Update auxiliary field values to current time.
     *
     * @param[in] t Current time.
     */
    void updateState(const PylithReal t);

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

    // PRIVATE MEMBERS /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    std::vector<ResidualKernels> _kernelsRHSResidual; ///< kernels for RHS residual.
    std::vector<ResidualKernels> _kernelsLHSResidual; ///< kernels for LHS residual.

    std::vector<JacobianKernels> _kernelsLHSJacobian; /// > kernels for LHS Jacobian.

    pylith::topology::Mesh* _interfaceMesh; ///< Boundary mesh.
    std::string _interfaceSurfaceLabel; ///< Name of label identifying interface surface.

    // NOT IMPLEMENTED /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    IntegratorInterface(void); ///< Not implemented.
    IntegratorInterface(const IntegratorInterface&); ///< Not implemented.
    const IntegratorInterface& operator=(const IntegratorInterface&); ///< Not implemented.

}; // IntegratorInterface

#endif // pylith_feassemble_integratorinterface_hh

// End of file
