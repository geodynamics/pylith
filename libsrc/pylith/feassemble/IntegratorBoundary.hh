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
 * @file libsrc/feassemble/IntegratorBoundary.hh
 *
 * @brief Object for finite-element integration over a boundary of the simulation domain.
 */

#if !defined(pylith_feassemble_integratorboundary_hh)
#define pylith_feassemble_integratorboundary_hh

#include "feassemblefwd.hh" // forward declarations

#include "pylith/feassemble/Integrator.hh" // ISA Integrator

#include "pylith/utils/arrayfwd.hh" // HASA std::vector

class pylith::feassemble::IntegratorBoundary : public pylith::feassemble::Integrator {
    friend class TestIntegratorBoundary; // unit testing

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

    // PUBLIC MEMBERS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    /// Constructor
    IntegratorBoundary(pylith::problems::Physics* const physics);

    /// Destructor
    virtual ~IntegratorBoundary(void);

    /// Deallocate PETSc and local data structures.
    virtual
    void deallocate(void);

    /** Set label marking boundary associated with boundary condition surface.
     *
     * @param[in] value Label of surface (from mesh generator).
     */
    void setMarkerLabel(const char* value);

    /** Get label marking boundary associated with boundary condition surface.
     *
     * @returns Label of surface (from mesh generator).
     */
    const char* getMarkerLabel(void) const;

    /** Get mesh associated with integrator domain.
     *
     * @returns Mesh associated with integrator domain.
     */
    const pylith::topology::Mesh& getIntegrationDomainMesh(void) const;

    /** Set kernels for RHS residual.
     *
     * @param kernels Array of kernerls for computing the RHS residual.
     */
    void setKernelsRHSResidual(const std::vector<ResidualKernels>& kernels);

    /** Set kernels for LHS residual.
     *
     * @param kernels Array of kernerls for computing the LHS residual.
     */
    void setKernelsLHSResidual(const std::vector<ResidualKernels>& kernels);

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
    std::vector<ResidualKernels> _kernelsLHSResidual; ///< kernels for LHS residual.

    pylith::topology::Mesh* _boundaryMesh; ///< Boundary mesh.
    std::string _boundaryLabel; ///< Label to identify boundary condition points in mesh.

    // NOT IMPLEMENTED /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    IntegratorBoundary(void); ///< Not implemented.
    IntegratorBoundary(const IntegratorBoundary&); ///< Not implemented.
    const IntegratorBoundary& operator=(const IntegratorBoundary&); ///< Not implemented.

}; // IntegratorBoundary

#endif // pylith_feassemble_integratorboundary_hh

// End of file
