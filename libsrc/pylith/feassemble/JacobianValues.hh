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
 * @file libsrc/feassemble/JacobianValues.hh
 *
 * @brief Object for setting Jacobian values without finite-element integration for a material.
 */

#if !defined(pylith_feassemble_jacobianvalues_hh)
#define pylith_feassemble_jacobianvalues_hh

#include "feassemblefwd.hh" // forward declarations
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
            function(NULL)
        {}


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

#endif // pylith_feassemble_jacobianvalues_hh

// End of file
