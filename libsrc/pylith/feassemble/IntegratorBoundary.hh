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
        EquationPart part; ///< Residual part (LHS or RHS).
        PetscBdPointFunc r0; ///< f0 (RHS) or g0 (LHS) function.
        PetscBdPointFunc r1; ///< f1 (RHS) or g1 (LHS) function.

        ResidualKernels(void) :
            subfield(""),
            part(pylith::feassemble::Integrator::LHS),
            r0(NULL),
            r1(NULL) {}


        ResidualKernels(const char* subfieldValue,
                        const EquationPart partValue,
                        PetscBdPointFunc r0Value,
                        PetscBdPointFunc r1Value) :
            subfield(subfieldValue),
            part(partValue),
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

    // PRIVATE MEMBERS /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    pylith::topology::Mesh* _boundaryMesh; ///< Boundary mesh.
    std::string _boundarySurfaceLabel; ///< Name of label identifying boundary surface.
    std::string _subfieldName; ///< Name of solution subfield for boundary condition.

    // NOT IMPLEMENTED /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    IntegratorBoundary(void); ///< Not implemented.
    IntegratorBoundary(const IntegratorBoundary&); ///< Not implemented.
    const IntegratorBoundary& operator=(const IntegratorBoundary&); ///< Not implemented.

}; // IntegratorBoundary

#endif // pylith_feassemble_integratorboundary_hh

// End of file
