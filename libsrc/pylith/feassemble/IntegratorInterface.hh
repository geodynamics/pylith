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
#include "pylith/feassemble/FEKernelKey.hh" // HASA FEKernelKey
#include "pylith/materials/materialsfwd.hh" // USES Material
#include "pylith/utils/arrayfwd.hh" // HASA std::vector

class pylith::feassemble::IntegratorInterface : public pylith::feassemble::Integrator {
    friend class _IntegratorInterface; // private utility class
    friend class TestIntegratorInterface; // unit testing

public:

    enum FaceEnum {
        NEGATIVE_FACE=0,
        POSITIVE_FACE=1,
        COHESIVE_FACE=2,
    }; // FaceEnum

    // PUBLIC STRUCTS /////////////////////////////////////////////////////////////////////////////
public:

    /// Kernels (pointwise functions) for residual.
    struct ResidualKernels {
        std::string subfield; ///< Name of subfield
        ResidualPart part; ///< Residual part (LHS or RHS).
        FaceEnum face; ///< Face domain.
        PetscBdPointFunc r0; ///< f0 (RHS) or g0 (LHS) function.
        PetscBdPointFunc r1; ///< f1 (RHS) or g1 (LHS) function.

        ResidualKernels(void) :
            subfield(""),
            part(pylith::feassemble::Integrator::RESIDUAL_LHS),
            face(COHESIVE_FACE),
            r0(NULL),
            r1(NULL) {}


        ResidualKernels(const char* subfieldValue,
                        const ResidualPart partValue,
                        FaceEnum faceValue,
                        PetscBdPointFunc r0Value,
                        PetscBdPointFunc r1Value) :
            subfield(subfieldValue),
            part(partValue),
            face(faceValue),
            r0(r0Value),
            r1(r1Value) {}


    }; // ResidualKernels

    /// Kernels (point-wise functions) for Jacobian;
    struct JacobianKernels {
        std::string subfieldTrial; ///< Name of subfield associated with trial function (row in Jacobian).
        std::string subfieldBasis; ///< Name of subfield associated with basis function (column in Jacobian).
        JacobianPart part; ///< Jacobian part (LHS or LHS lumped inverse).
        FaceEnum face; ///< Integration domain.
        PetscBdPointJac j0; ///< J0 function.
        PetscBdPointJac j1; ///< J1 function.
        PetscBdPointJac j2; ///< J2 function.
        PetscBdPointJac j3; ///< J3 function.

        JacobianKernels(void) :
            subfieldTrial(""),
            subfieldBasis(""),
            part(JACOBIAN_LHS),
            face(COHESIVE_FACE),
            j0(NULL),
            j1(NULL),
            j2(NULL),
            j3(NULL) {}


        JacobianKernels(const char* subfieldTrialValue,
                        const char* subfieldBasisValue,
                        JacobianPart partValue,
                        FaceEnum faceValue,
                        PetscBdPointJac j0Value,
                        PetscBdPointJac j1Value,
                        PetscBdPointJac j2Value,
                        PetscBdPointJac j3Value) :
            subfieldTrial(subfieldTrialValue),
            subfieldBasis(subfieldBasisValue),
            part(partValue),
            face(faceValue),
            j0(j0Value),
            j1(j1Value),
            j2(j2Value),
            j3(j3Value) {}


    }; // JacobianKernels

    // PUBLIC MEMBERS /////////////////////////////////////////////////////////////////////////////
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

    /** Set integration patches.
     *
     * @param[in] patches Interface integration patches.
     */
    void setIntegrationPatches(pylith::feassemble::InterfacePatches* patches);

    /** Set kernels for residual.
     *
     * @param kernels Array of kernels for computing the residual.
     * @param[in] solution Field with current trial solution.
     */
    void setKernels(const std::vector<ResidualKernels>& kernels,
                    const pylith::topology::Field& solution,
                    const std::vector<pylith::materials::Material*>& materials);

    /** Set kernels for Jacobian.
     *
     * @param kernels Array of kernels for computing the Jacobian.
     * @param[in] solution Field with current trial solution.
     */
    void setKernels(const std::vector<JacobianKernels>& kernels,
                    const pylith::topology::Field& solution,
                    const std::vector<pylith::materials::Material*>& materials);

    /** Initialize integration domain, auxiliary field, and derived field. Update observers.
     *
     * @param[in] solution Solution field (layout).
     */
    void initialize(const pylith::topology::Field& solution);

    /** Set auxiliary field values for current time.
     *
     * @param[in] t Current time.
     * @param[in] dt Current time step.
     */
    void setState(const PylithReal t,
                  const PylithReal dt);

    /** Compute RHS residual for G(t,s).
     *
     * @param[out] residual Field for residual.
     * @param[in] integrationData Data needed to integrate governing equations.
     */
    void computeRHSResidual(pylith::topology::Field* residual,
                            const pylith::problems::IntegrationData& integrationData);

    /** Compute LHS residual for F(t,s,\dot{s}).
     *
     * @param[out] residual Field for residual.
     * @param[in] integrationData Data needed to integrate governing equations.
     */
    void computeLHSResidual(pylith::topology::Field* residual,
                            const pylith::problems::IntegrationData& integrationData);

    /** Compute LHS Jacobian and preconditioner for F(t,s,\dot{s}) with implicit time-stepping.
     *
     * @param[out] jacobianMat PETSc Mat with Jacobian sparse matrix.
     * @param[out] precondMat PETSc Mat with Jacobian preconditioning sparse matrix.
     * @param[in] integrationData Data needed to integrate governing equations.
     */
    void computeLHSJacobian(PetscMat jacobianMat,
                            PetscMat precondMat,
                            const pylith::problems::IntegrationData& integrationData);

    /** Compute inverse of lumped LHS Jacobian for F(t,s,\dot{s}) with explicit time-stepping.
     *
     * @param[out] jacobianInv Inverse of lumped Jacobian as a field.
     * @param[in] integrationData Data needed to integrate governing equations.
     */
    void computeLHSJacobianLumpedInv(pylith::topology::Field* jacobianInv,
                                     const pylith::problems::IntegrationData& integrationData);

    // PRIVATE ENUM ///////////////////////////////////////////////////////////////////////////////
private:

    enum ResidualPartKey {
        RESIDUAL_LHS_NEG=0, // LHS residual negative face.
        RESIDUAL_LHS_POS=1, // LHS residual positive face.
        RESIDUAL_LHS_FAULT=2, // LHS residual hybrid cell.
        RESIDUAL_RHS_NEG=4, // RHS residual negative face.
        RESIDUAL_RHS_POS=5, // RHS residual negative face.
        RESIDUAL_RHS_FAULT=6, // RHS residual negative face.
        RESIDUAL_LHS_WT_NEG=7, // LHS weighted residual negative face.
        RESIDUAL_LHS_WT_POS=8, // LHS weighted residual negative face.
        RESIDUAL_LHS_WT_FAULT=9, // LHS weighted residual negative face.
    };

    enum JacobianPartKey {
        JACOBIAN_LHS_NEG=0, // LHS Jacobian negative face.
        JACOBIAN_LHS_POS=1, // LHS Jacobian positive face.
        JACOBIAN_LHS_FAULT=2, // LHS Jacobian hybrid cell.
    };

    // PRIVATE MEMBERS ////////////////////////////////////////////////////////////////////////////
private:

    pylith::topology::Mesh* _interfaceMesh; ///< Boundary mesh.
    std::string _interfaceSurfaceLabel; ///< Name of label identifying interface surface.

    pylith::feassemble::InterfacePatches* _integrationPatches; ///< Face patches.

    PetscDM _weightingDM; ///< PETSc DM for weighting.
    PetscVec _weightingVec; ///< PETSc Vec for weighting values.

    // NOT IMPLEMENTED ////////////////////////////////////////////////////////////////////////////
private:

    IntegratorInterface(void); ///< Not implemented.
    IntegratorInterface(const IntegratorInterface&); ///< Not implemented.
    const IntegratorInterface& operator=(const IntegratorInterface&); ///< Not implemented.

}; // IntegratorInterface

#endif // pylith_feassemble_integratorinterface_hh

// End of file
