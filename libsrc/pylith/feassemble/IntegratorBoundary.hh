// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================
#pragma once

#include "pylith/feassemble/feassemblefwd.hh" // forward declarations

#include "pylith/feassemble/Integrator.hh" // ISA Integrator

#include "pylith/utils/arrayfwd.hh" // HASA std::vector

class pylith::feassemble::IntegratorBoundary : public pylith::feassemble::Integrator {
    friend class TestIntegratorBoundary; // unit testing

    // PUBLIC STRUCTS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    /// Kernels (point-wise functions) for residual.
    struct ResidualKernels {
        std::string subfield; ///< Name of subfield
        EquationPart part; ///< Residual part (LHS or RHS).
        PetscBdPointFn* r0; ///< f0 (RHS) or g0 (LHS) function.
        PetscBdPointFn* r1; ///< f1 (RHS) or g1 (LHS) function.

        ResidualKernels(void) :
            subfield(""),
            part(pylith::feassemble::Integrator::LHS),
            r0(NULL),
            r1(NULL) {}


        ResidualKernels(const char* subfieldValue,
                        const EquationPart partValue,
                        PetscBdPointFn* r0Value,
                        PetscBdPointFn* r1Value) :
            subfield(subfieldValue),
            part(partValue),
            r0(r0Value),
            r1(r1Value) {}


    }; // ResidualKernels

    /// Project kernels (pointwise functions) for updating state variables or computing derived fields.
    struct ProjectKernels {
        std::string subfield; ///< Name of subfield for function.
        PetscBdPointFn* f; ///< Point-wise function.

        ProjectKernels(void) :
            subfield(""),
            f(NULL) {}


        ProjectKernels(const char* subfieldValue,
                       PetscBdPointFn* fValue) :
            subfield(subfieldValue),
            f(fValue) {}


    }; // ProjectKernels

    // PUBLIC MEMBERS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    /// Constructor
    IntegratorBoundary(pylith::problems::Physics* const physics);

    /// Destructor
    virtual ~IntegratorBoundary(void);

    /// Deallocate PETSc and local data structures.
    virtual
    void deallocate(void);

    /** Set name of solution subfield associated with boundary condition.
     *
     * @param[in] value Name of solution subfield.
     */
    void setSubfieldName(const char* value);

    /** Get name of solution subfield associated with boundary condition.
     *
     * @preturn Name of solution subfield.
     */
    const char* getSubfieldName(void) const;

    /** Get mesh associated with integrator domain.
     *
     * @returns Mesh associated with integrator domain.
     */
    const pylith::topology::Mesh& getPhysicsDomainMesh(void) const;

    /** Set kernels for RHS residual.
     *
     * @param[in] kernels Array of kernerls for computing the RHS residual.
     * @param[in] solution Solution field.
     */
    void setKernelsResidual(const std::vector<ResidualKernels>& kernels,
                            const pylith::topology::Field& solution);

    /** Set kernels for computing diagnostic field.
     *
     * @param kernels Array of kernels for computing diagnostic field.
     */
    void setKernelsDiagnosticField(const std::vector<ProjectKernels>& kernels);

    /** Initialize integration domain, auxiliary field, and derived field. Update observers.
     *
     * @param[in] solution Solution field (layout).
     */
    void initialize(const pylith::topology::Field& solution);

    /** Set auxiliary field for current time.
     *
     * @param[in] t Current time.
     */
    void setState(const PylithReal t);

    /** Compute RHS residual for G(t,s).
     *
     * @param[out] residual Field for residual.
     * @param[in] integrationData Data needed to integrate governing equations.
     */
    void computeRHSResidual(pylith::topology::Field* residual,
                            const pylith::feassemble::IntegrationData& integrationData);

    /** Compute LHS residual for F(t,s,\dot{s}).
     *
     * @param[out] residual Field for residual.
     * @param[in] integrationData Data needed to integrate governing equations.
     */
    void computeLHSResidual(pylith::topology::Field* residual,
                            const pylith::feassemble::IntegrationData& integrationData);

    /** Compute LHS Jacobian and preconditioner for F(t,s,\dot{s}) with implicit time-stepping.
     *
     * @param[out] jacobianMat PETSc Mat with Jacobian sparse matrix.
     * @param[out] precondMat PETSc Mat with Jacobian preconditioning sparse matrix.
     * @param[in] integrationData Data needed to integrate governing equations.
     */
    void computeLHSJacobian(PetscMat jacobianMat,
                            PetscMat precondMat,
                            const pylith::feassemble::IntegrationData& integrationData);

    /** Compute inverse of lumped LHS Jacobian for F(t,s,\dot{s}) with explicit time-stepping.
     *
     * @param[out] jacobianInv Inverse of lumped Jacobian as a field.
     * @param[in] integrationData Data needed to integrate governing equations.
     */
    void computeLHSJacobianLumpedInv(pylith::topology::Field* jacobianInv,
                                     const pylith::feassemble::IntegrationData& integrationData);

    // PROTECTED METHODS //////////////////////////////////////////////////////////////////////////
protected:

    /// Compute diagnostic field from auxiliary field.
    virtual
    void _computeDiagnosticField(void);

    // PRIVATE MEMBERS /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    pylith::topology::Mesh* _boundaryMesh; ///< Boundary mesh.
    std::string _boundarySurfaceLabel; ///< Name of label identifying boundary surface.
    std::string _subfieldName; ///< Name of solution subfield for boundary condition.

    std::vector<ProjectKernels> _kernelsDiagnosticField; ///< kernels for computing diagnostic field.

    // NOT IMPLEMENTED /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    IntegratorBoundary(void); ///< Not implemented.
    IntegratorBoundary(const IntegratorBoundary&); ///< Not implemented.
    const IntegratorBoundary& operator=(const IntegratorBoundary&); ///< Not implemented.

}; // IntegratorBoundary

// End of file
