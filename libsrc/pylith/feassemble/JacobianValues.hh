// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2024, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================
#pragma once

#include "pylith/feassemble/feassemblefwd.hh" // forward declarations
#include "pylith/utils/GenericComponent.hh" // ISA GenericComponent

#include "pylith/topology/topologyfwd.hh" // USES Field
#include "pylith/utils/petscfwd.h" // USES PetscMat
#include "pylith/utils/arrayfwd.hh" // HASA std::vector

class pylith::feassemble::JacobianValues : public pylith::utils::GenericComponent {
    friend class TestJacobianValues; // unit testing

    // PUBLIC STRUCTS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    typedef void (*PointFn)(pylith::scalar_array*,
                            const PylithReal,
                            const PylithReal,
                            const PylithReal,
                            const PylithInt,
                            const PylithInt,
                            const PylithInt,
                            const PylithInt,
                            const PylithInt);

    /// Kernels (point-wise functions) for Jacobian;
    struct JacobianKernel {
        std::string subfieldTrial; ///< Name of subfield associated with trial function (row in Jacobian).
        std::string subfieldBasis; ///< Name of subfield associated with basis function (column in Jacobian).
        PointFn function; ///< Pointwise function for Jacobian values.

        JacobianKernel(void) :
            subfieldTrial(""),
            subfieldBasis(""),
            function(NULL) {}


        JacobianKernel(const char* subfieldTrialValue,
                       const char* subfieldBasisValue,
                       const PointFn functionValue) :
            subfieldTrial(subfieldTrialValue),
            subfieldBasis(subfieldBasisValue),
            function(functionValue) {}


    }; // JacobianKernel

    // PUBLIC MEMBERS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    /// Constructor
    JacobianValues(void);

    /// Destructor
    ~JacobianValues(void);

    /** Set kernels.
     *
     * @param[in] kernelsJacobian Array of kernels for computing the Jacobian.
     * @param[in] kernelsPrecond Array of kernels for computing the preconditioner.
     */
    void setKernels(const std::vector<JacobianKernel>& kernelsJacobian,
                    const std::vector<JacobianKernel>& kernelsPrecond);

    /** Compute LHS Jacobian and preconditioner for F(t,s,\dot{s}) with implicit time-stepping.
     *
     * @param[out] jacobianMat PETSc Mat with Jacobian sparse matrix.
     * @param[out] precondMat PETSc Mat with Jacobian preconditioning sparse matrix.
     * @param[in] t Current time.
     * @param[in] dt Current time step.
     * @param[in] s_tshift Scale for time derivative.
     * @param[in] solution Field with current trial solution.
     * @param[in] dsLabel PETSc DS label information.
     */
    void computeLHSJacobian(PetscMat jacobianMat,
                            PetscMat precondMat,
                            const PylithReal t,
                            const PylithReal dt,
                            const PylithReal s_tshift,
                            const pylith::topology::Field& solution,
                            const pylith::feassemble::DSLabelAccess& dsLabel);

    /** Compute Jacobian values associated with s_tshift on diagonal.
     *
     * @param[out] jacobianMat PETSc Mat with Jacobian sparse matrix.
     * @param[in] t Current time.
     * @param[in] dt Current time step.
     * @param[in] s_tshift Scale for time derivative.
     */
    static
    void blockDiag_tshift(pylith::scalar_array* cellMat,
                          const PylithReal t,
                          const PylithReal dt,
                          const PylithReal s_tshift,
                          const PylithInt trialDof,
                          const PylithInt trialOff,
                          const PylithInt basisDof,
                          const PylithInt basisOff,
                          const PylithInt totalDim);

    // PRIVATE MEMBERS /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    std::vector<JacobianKernel> _kernelsJacobian; ///< Kernels for Jacobian.
    std::vector<JacobianKernel> _kernelsPrecond; ///< Kernels for preconditioner.

    // NOT IMPLEMENTED /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    JacobianValues(const JacobianValues&); ///< Not implemented.
    const JacobianValues& operator=(const JacobianValues&); ///< Not implemented.

}; // JacobianValues

// End of file
