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
// Copyright (c) 2010-2022 University of California, Davis
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
        FAULT_FACE=2,
    }; // FaceEnum

    // PUBLIC STRUCTS /////////////////////////////////////////////////////////////////////////////
public:

    /// Kernels (pointwise functions) for residual.
    struct ResidualKernels {
        std::string subfield; ///< Name of subfield
        EquationPart part; ///< Residual part (LHS or RHS).
        FaceEnum face; ///< Face domain.
        PetscBdPointFunc r0; ///< f0 (RHS) or g0 (LHS) function.
        PetscBdPointFunc r1; ///< f1 (RHS) or g1 (LHS) function.

        ResidualKernels(void) :
            subfield(""),
            part(pylith::feassemble::Integrator::LHS),
            face(FAULT_FACE),
            r0(NULL),
            r1(NULL) {}


        ResidualKernels(const char* subfieldValue,
                        const EquationPart partValue,
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
        EquationPart part; ///< Jacobian part (LHS or LHS lumped inverse).
        FaceEnum face; ///< Integration domain.
        PetscBdPointJac j0; ///< J0 function.
        PetscBdPointJac j1; ///< J1 function.
        PetscBdPointJac j2; ///< J2 function.
        PetscBdPointJac j3; ///< J3 function.

        JacobianKernels(void) :
            subfieldTrial(""),
            subfieldBasis(""),
            part(LHS),
            face(FAULT_FACE),
            j0(NULL),
            j1(NULL),
            j2(NULL),
            j3(NULL) {}


        JacobianKernels(const char* subfieldTrialValue,
                        const char* subfieldBasisValue,
                        EquationPart partValue,
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

    /** Set name of label marking interior interface surface.
     *
     * @param[in] value Name of label of surface (from mesh generator).
     */
    void setSurfaceLabelName(const char* value);

    /** Get name of label marking interior interface surface.
     *
     * @returns Name of label of surface (from mesh generator).
     */
    const char* getSurfaceLabelName(void) const;

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

    /** Get integration patches.
     *
     * @returns Interface integration patches.
     */
    const pylith::feassemble::InterfacePatches* getIntegrationPatches(void) const;

    /** Set kernels for residual.
     *
     * @param kernels Array of kernels for computing the residual.
     * @param[in] solution Field with current trial solution.
     * @param[in] materials Materials in problem.
     */
    void setKernels(const std::vector<ResidualKernels>& kernels,
                    const pylith::topology::Field& solution,
                    const std::vector<pylith::materials::Material*>& materials);

    /** Set kernels for Jacobian.
     *
     * @param kernels Array of kernels for computing the Jacobian.
     * @param[in] solution Field with current trial solution.
     * @param[in] materials Materials in problem.
     */
    void setKernels(const std::vector<JacobianKernels>& kernels,
                    const pylith::topology::Field& solution,
                    const std::vector<pylith::materials::Material*>& materials);

    /** Compute weak form key part for face.
     *
     * For integration with hybrid cells, we must distinguish among integration of the
     * negative and positive faces for each fault, the fault face as well as the residual term.
     *
     * @param[in] eqnPart Term in the equation.
     * @param[in] face Negative, positive, or fault face.
     * @param[in] patch Interface patch label value.
     */
    PetscInt getWeakFormPart(const PetscInt eqnPart,
                             const PetscInt face,
                             const PetscInt patch) const;

    /** Initialize integration domain, auxiliary field, and derived field. Update observers.
     *
     * @param[in] solution Solution field (layout).
     */
    void initialize(const pylith::topology::Field& solution);

    /** Set auxiliary field values for current time.
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

    // PRIVATE MEMBERS ////////////////////////////////////////////////////////////////////////////
private:

    pylith::topology::Mesh* _interfaceMesh; ///< Boundary mesh.
    std::string _surfaceLabelName; ///< Name of label identifying interface surface.

    pylith::feassemble::InterfacePatches* _integrationPatches; ///< Face patches.

    PetscDM _weightingDM; ///< PETSc DM for weighting.
    PetscVec _weightingVec; ///< PETSc Vec for weighting values.

    bool _hasLHSResidualWeighted; ///< Has LHS Residual with weighted terms.
    bool _hasLHSJacobianWeighted; ///< Has LHS Jacobian with weighted terms.

    // NOT IMPLEMENTED ////////////////////////////////////////////////////////////////////////////
private:

    IntegratorInterface(void); ///< Not implemented.
    IntegratorInterface(const IntegratorInterface&); ///< Not implemented.
    const IntegratorInterface& operator=(const IntegratorInterface&); ///< Not implemented.

}; // IntegratorInterface

#endif // pylith_feassemble_integratorinterface_hh

// End of file
