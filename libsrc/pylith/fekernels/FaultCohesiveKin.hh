/*
 * ================================================================================================
 * This code is part of PyLith, developed through the Computational Infrastructure
 * for Geodynamics (https://github.com/geodynamics/pylith).
 *
 * Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
 * All rights reserved.
 *
 * See https://mit-license.org/ and LICENSE.md and for license information.
 * =================================================================================================
 */
#pragma once

/*
 * Kernels for faults with prescribed slip.
 *
 * Solution fields: [disp(dim), vel(dim, optional), lagrange(dim)]
 *
 * Auxiliary fields: [...] (not used)
 *
 * LHS Residual
 *
 *  - f0u^+ = -\lambda
 *  - f0u^- = +\lambda
 *  - f0\lambda = d - u^+ + u^-
 *
 * LHS Jacobian
 *
 *  - Jf0^{u \lambda} = -1 for u^+, +1 for u^-
 *  - Jf0^{\lambda u} = -1 for u^+, +1 for u^-
 *
 * Kernel interface.
 *
 * @param[in] dim Spatial dimension.
 * @param[in] numS Number of registered subfields in solution field.
 * @param[in] numA Number of registered subfields in auxiliary field.
 * @param[in] sOff Offset of registered subfields in solution field [numS].
 * @param[in] sOff_x Offset of registered subfields in gradient of the solution field [numS].
 * @param[in] s Solution field with all subfields.
 * @param[in] s_t Time derivative of solution field.
 * @param[in] s_x Gradient of solution field.
 * @param[in] aOff Offset of registered subfields in auxiliary field [numA]
 * @param[in] aOff_x Offset of registered subfields in gradient of auxiliary field [numA]
 * @param[in] a Auxiliary field with all subfields.
 * @param[in] a_t Time derivative of auxiliary field.
 * @param[in] a_x Gradient of auxiliary field.
 * @param[in] t Time for residual evaluation.
 * @param[in] x Coordinates of point evaluation.
 * @param[in] n Outward normal vector for boundary.
 * @param[in] numConstants Number of registered constants.
 * @param[in] constants Array of registered constants.
 * @param[out] f0 [dim].
 */

#include "pylith/fekernels/fekernelsfwd.hh" // forward declarations

#include "pylith/fekernels/Tensor.hh" // USES Tensor
#include "pylith/fekernels/BoundaryDirections.hh" // USES BoundaryDirections
#include "pylith/fekernels/Elasticity.hh" // USES Elasticity

#include "pylith/utils/types.hh"
#include "pylith/utils/journals.hh"

class pylith::fekernels::FaultCohesiveKin {
    // PUBLIC METHODS /////////////////////////////////////////////////////////////////////////////
public:

    // --------------------------------------------------------------------------------------------
    /** f0 function for elasticity equation: f0u = +\lambda (neg side).
     *
     * Solution fields: [disp(dim), ..., lagrange(dim)]
     */
    static inline
    void f0u_neg(const PylithInt dim,
                 const PylithInt numS,
                 const PylithInt numA,
                 const PylithInt sOff[],
                 const PylithInt sOff_x[],
                 const PylithScalar s[],
                 const PylithScalar s_t[],
                 const PylithScalar s_x[],
                 const PylithInt aOff[],
                 const PylithInt aOff_x[],
                 const PylithScalar a[],
                 const PylithScalar a_t[],
                 const PylithScalar a_x[],
                 const PylithReal t,
                 const PylithScalar x[],
                 const PylithReal n[],
                 const PylithInt numConstants,
                 const PylithScalar constants[],
                 PylithScalar f0[]) {
        assert(sOff);
        assert(s);
        assert(f0);

        assert(numS >= 2);

        const PylithInt spaceDim = dim; // dim passed is spatial dimension

        const PylithInt fOffN = 0;
        const PylithInt sOffLagrange = sOff[numS-1];
        const PylithScalar* lagrange = &s[sOffLagrange];

        for (PylithInt i = 0; i < spaceDim; ++i) {
            f0[fOffN+i] += +lagrange[i];
        } // for
    }

    // --------------------------------------------------------------------------------------------
    /** f0 function for elasticity equation: f0u = -\lambda (pos side).
     *
     * Solution fields: [disp(dim), ..., lagrange(dim)]
     */
    static inline
    void f0u_pos(const PylithInt dim,
                 const PylithInt numS,
                 const PylithInt numA,
                 const PylithInt sOff[],
                 const PylithInt sOff_x[],
                 const PylithScalar s[],
                 const PylithScalar s_t[],
                 const PylithScalar s_x[],
                 const PylithInt aOff[],
                 const PylithInt aOff_x[],
                 const PylithScalar a[],
                 const PylithScalar a_t[],
                 const PylithScalar a_x[],
                 const PylithReal t,
                 const PylithScalar x[],
                 const PylithReal n[],
                 const PylithInt numConstants,
                 const PylithScalar constants[],
                 PylithScalar f0[]) {
        assert(sOff);
        assert(s);
        assert(f0);

        assert(numS >= 2);

        const PylithInt spaceDim = dim; // dim passed is spatial dimension

        const PylithInt fOffP = 0;
        const PylithInt sOffLagrange = sOff[numS-1];
        const PylithScalar* lagrange = &s[sOffLagrange];

        for (PylithInt i = 0; i < spaceDim; ++i) {
            f0[fOffP+i] += -lagrange[i];
        } // for
    }

    // --------------------------------------------------------------------------------------------
    /** f0 function for slip constraint equation: f0\lambda = (u^+ - u^-) - d
     *
     * Solution fields: [disp(dim), ..., lagrange(dim)]
     */
    static inline
    void f0l_slip(const PylithInt dim,
                  const PylithInt numS,
                  const PylithInt numA,
                  const PylithInt sOff[],
                  const PylithInt sOff_x[],
                  const PylithScalar s[],
                  const PylithScalar s_t[],
                  const PylithScalar s_x[],
                  const PylithInt aOff[],
                  const PylithInt aOff_x[],
                  const PylithScalar a[],
                  const PylithScalar a_t[],
                  const PylithScalar a_x[],
                  const PylithReal t,
                  const PylithScalar x[],
                  const PylithReal n[],
                  const PylithInt numConstants,
                  const PylithScalar constants[],
                  PylithScalar f0[]) {
        assert(sOff);
        assert(aOff);
        assert(s);
        assert(a);
        assert(f0);

        assert(numS >= 2);
        assert(numA >= 1);

        const PylithInt spaceDim = dim; // dim passed is spatial dimension
        const PylithInt i_slip = 0;
        const PylithInt i_disp = 0;

        const PylithScalar* slip = &a[aOff[i_slip]];

        const PylithInt sOffDispN = sOff[i_disp];
        const PylithInt sOffDispP = sOffDispN+spaceDim;
        const PylithInt fOffLagrange = 0;

        const PylithScalar* dispN = &s[sOffDispN];
        const PylithScalar* dispP = &s[sOffDispP];

        switch (spaceDim) {
        case 2: {
            const PylithInt _spaceDim = 2;
            const PylithScalar tanDir[2] = {-n[1], n[0] };
            for (PylithInt i = 0; i < _spaceDim; ++i) {
                const PylithScalar slipXY = n[i]*slip[0] + tanDir[i]*slip[1];
                f0[fOffLagrange+i] += -dispP[i] + dispN[i] + slipXY;
            } // for
            break;
        } // case 2
        case 3: {
            const PylithInt _spaceDim = 3;
            const PylithScalar* refDir1 = &constants[0];
            const PylithScalar* refDir2 = &constants[3];
            PylithScalar tanDir1[3], tanDir2[3];
            pylith::fekernels::BoundaryDirections::tangential_directions(tanDir1, tanDir2, refDir1, refDir2, n);

            for (PylithInt i = 0; i < _spaceDim; ++i) {
                const PylithScalar slipXYZ = n[i]*slip[0] + tanDir1[i]*slip[1] + tanDir2[i]*slip[2];
                f0[fOffLagrange+i] += -dispP[i] + dispN[i] + slipXYZ;
            } // for
            break;
        } // case 3
        default:
            assert(0);
        } // switch
    }

    // --------------------------------------------------------------------------------------------
    /** f0 function for slip acceleration constraint equation: f0\lambda = (\dot{v}^+ - \dot{v}^-) - \ddot{d}
     *
     * Solution fields: [disp(dim), vel(dim), ..., lagrange(dim)]
     */
    static inline
    void f0l_slipAcc(const PylithInt dim,
                     const PylithInt numS,
                     const PylithInt numA,
                     const PylithInt sOff[],
                     const PylithInt sOff_x[],
                     const PylithScalar s[],
                     const PylithScalar s_t[],
                     const PylithScalar s_x[],
                     const PylithInt aOff[],
                     const PylithInt aOff_x[],
                     const PylithScalar a[],
                     const PylithScalar a_t[],
                     const PylithScalar a_x[],
                     const PylithReal t,
                     const PylithScalar x[],
                     const PylithReal n[],
                     const PylithInt numConstants,
                     const PylithScalar constants[],
                     PylithScalar f0[]) {
        assert(sOff);
        assert(aOff);
        assert(s);
        assert(a);
        assert(f0);

        assert(numS >= 3);
        assert(numA >= 1);

        const PylithInt spaceDim = dim; // dim passed is spatial dimension
        const PylithInt i_slipAcc = numA-1;

        const PylithScalar* slipAcc = &a[aOff[i_slipAcc]];

        switch (spaceDim) {
        case 2: {
            const PylithInt _spaceDim = 2;
            const PylithScalar tanDir[2] = {-n[1], n[0] };
            for (PylithInt i = 0; i < _spaceDim; ++i) {
                const PylithScalar slipAccXY = n[i]*slipAcc[0] + tanDir[i]*slipAcc[1];
                f0[i] -= slipAccXY;
            } // for
            break;
        } // case 2
        case 3: {
            const PylithInt _spaceDim = 3;
            const PylithScalar* refDir1 = &constants[0];
            const PylithScalar* refDir2 = &constants[3];
            PylithScalar tanDir1[3], tanDir2[3];
            pylith::fekernels::BoundaryDirections::BoundaryDirections::tangential_directions(tanDir1, tanDir2, refDir1, refDir2, n);

            for (PylithInt i = 0; i < _spaceDim; ++i) {
                const PylithScalar slipAccXYZ = n[i]*slipAcc[0] + tanDir1[i]*slipAcc[1] + tanDir2[i]*slipAcc[2];
                f0[i] -= slipAccXYZ;
            } // for
            break;
        } // case 3
        default:
            assert(0);
        } // switch
    }

    // --------------------------------------------------------------------------------------------
    /** Jf0 function for displacement equation: -\lambda (neg side).
     */
    static inline
    void Jf0ul_neg(const PylithInt dim,
                   const PylithInt numS,
                   const PylithInt numA,
                   const PylithInt sOff[],
                   const PylithInt sOff_x[],
                   const PylithScalar s[],
                   const PylithScalar s_t[],
                   const PylithScalar s_x[],
                   const PylithInt aOff[],
                   const PylithInt aOff_x[],
                   const PylithScalar a[],
                   const PylithScalar a_t[],
                   const PylithScalar a_x[],
                   const PylithReal t,
                   const PylithReal s_tshift,
                   const PylithScalar x[],
                   const PylithReal n[],
                   const PylithInt numConstants,
                   const PylithScalar constants[],
                   PylithScalar Jf0[]) {
        assert(numS >= 2);
        assert(Jf0);
        assert(sOff);
        assert(n);

        const PylithInt spaceDim = dim; // dim passed is spatial dimension

        const PylithInt gOffN = 0;
        const PylithInt ncols = spaceDim;

        for (PylithInt i = 0; i < spaceDim; ++i) {
            Jf0[(gOffN+i)*ncols+i] += +1.0;
        } // for
    }

    // --------------------------------------------------------------------------------------------
    /** Jf0 function for displacement equation: +\lambda (pos side).
     */
    static inline
    void Jf0ul_pos(const PylithInt dim,
                   const PylithInt numS,
                   const PylithInt numA,
                   const PylithInt sOff[],
                   const PylithInt sOff_x[],
                   const PylithScalar s[],
                   const PylithScalar s_t[],
                   const PylithScalar s_x[],
                   const PylithInt aOff[],
                   const PylithInt aOff_x[],
                   const PylithScalar a[],
                   const PylithScalar a_t[],
                   const PylithScalar a_x[],
                   const PylithReal t,
                   const PylithReal s_tshift,
                   const PylithScalar x[],
                   const PylithReal n[],
                   const PylithInt numConstants,
                   const PylithScalar constants[],
                   PylithScalar Jf0[]) {
        assert(numS >= 2);
        assert(Jf0);
        assert(sOff);
        assert(n);

        const PylithInt spaceDim = dim; // dim passed is spatial dimension

        const PylithInt ncols = spaceDim;

        for (PylithInt i = 0; i < spaceDim; ++i) {
            Jf0[i*ncols+i] += -1.0;
        } // for
    }

    // --------------------------------------------------------------------------------------------
    /** Jf0 function for slip constraint equation: +\lambda (pos side), -\lambda (neg side).
     *
     * Solution fields: [disp(dim), ..., lagrange(dim)]
     */
    static inline
    void Jf0lu(const PylithInt dim,
               const PylithInt numS,
               const PylithInt numA,
               const PylithInt sOff[],
               const PylithInt sOff_x[],
               const PylithScalar s[],
               const PylithScalar s_t[],
               const PylithScalar s_x[],
               const PylithInt aOff[],
               const PylithInt aOff_x[],
               const PylithScalar a[],
               const PylithScalar a_t[],
               const PylithScalar a_x[],
               const PylithReal t,
               const PylithReal s_tshift,
               const PylithScalar x[],
               const PylithReal n[],
               const PylithInt numConstants,
               const PylithScalar constants[],
               PylithScalar Jf0[]) {
        assert(numS >= 2);
        assert(Jf0);
        assert(sOff);
        assert(n);

        const PylithInt spaceDim = dim; // dim passed is spatial dimension

        const PylithInt gOffN = 0;
        const PylithInt gOffP = gOffN+spaceDim*spaceDim;
        const PylithInt ncols = spaceDim;

        for (PylithInt i = 0; i < spaceDim; ++i) {
            Jf0[gOffN+i*ncols+i] += +1.0;
            Jf0[gOffP+i*ncols+i] += -1.0;
        } // for
    }

    // ------------------------------------------------------------------------------------------------
    /** f0 function for negative fault face for DYNAMIC_IMEX formulation.
     *
     *  Solution fields: [disp(dim), vel(dim), lagrange(dim)]
     */
    static inline
    void f0l_neg(const PylithInt dim,
                 const PylithInt numS,
                 const PylithInt sOff[],
                 const PylithScalar s[],
                 const PylithReal n[],
                 const pylith::fekernels::Elasticity::StrainContext& strainContext,
                 void* rheologyContext,
                 pylith::fekernels::Elasticity::strainfn_type strainFn,
                 pylith::fekernels::Elasticity::stressfn_type stressFn,
                 pylith::fekernels::Elasticity::tractionfn_type tractionFn,
                 const pylith::fekernels::TensorOps& tensorOps,
                 PylithScalar f0[]) {
        const PylithInt spaceDim = dim; // dim passed is spatial dimension

        // Incoming solution fields.
        const PylithInt i_lagrange = 2;

        assert(numS >= 1);
        assert(sOff);
        assert(sOff[i_lagrange] >= 0);
        assert(f0);

        pylith::fekernels::Tensor strain;
        strainFn(strainContext, &strain);

        pylith::fekernels::Tensor stress;
        stressFn(rheologyContext, strain, tensorOps, &stress);

        PylithReal traction[3] = { 0.0, 0.0, 0.0 };
        tractionFn(stress, n, traction);

        const PylithScalar* lagrange = &s[sOff[i_lagrange]];
        for (PylithInt i = 0; i < spaceDim; ++i) {
            f0[i] += lagrange[i] - traction[i];
        } // for
    } // f0l_neg

    // ------------------------------------------------------------------------------------------------
    /** f0 function for positive fault face for DYNAMIC_IMEX formulation.
     *
     *  Solution fields: [disp(dim), vel(dim), lagrange(dim)]
     */
    static inline
    void f0l_pos(const PylithInt dim,
                 const PylithInt numS,
                 const PylithInt sOff[],
                 const PylithScalar s[],
                 const PylithReal n[],
                 const pylith::fekernels::Elasticity::StrainContext& strainContext,
                 void* rheologyContext,
                 pylith::fekernels::Elasticity::strainfn_type strainFn,
                 pylith::fekernels::Elasticity::stressfn_type stressFn,
                 pylith::fekernels::Elasticity::tractionfn_type tractionFn,
                 const pylith::fekernels::TensorOps& tensorOps,
                 PylithScalar f0[]) {
        const PylithInt spaceDim = dim; // dim passed is spatial dimension

        // Incoming solution fields.
        const PylithInt i_lagrange = 2;

        assert(numS >= 1);
        assert(sOff);
        assert(sOff[i_lagrange] >= 0);
        assert(f0);

        pylith::fekernels::Tensor strain;
        strainFn(strainContext, &strain);

        pylith::fekernels::Tensor stress;
        stressFn(rheologyContext, strain, tensorOps, &stress);

        PylithReal traction[3] = { 0.0, 0.0, 0.0 };
        tractionFn(stress, n, traction);

        const PylithScalar* lagrange = &s[sOff[i_lagrange]];
        for (PylithInt i = 0; i < spaceDim; ++i) {
            f0[i] += lagrange[i] + traction[i];
        } // for
    } // f0l_pos

    // ------------------------------------------------------------------------------------------------
    // Jf0 function for dynamic slip constraint equation for negative side of the fault.
    static inline
    void Jf0ll_neg(const PylithInt dim,
                   const PylithInt numS,
                   const PylithInt numA,
                   const PylithInt sOff[],
                   const PylithInt sOff_x[],
                   const PylithScalar s[],
                   const PylithScalar s_t[],
                   const PylithScalar s_x[],
                   const PylithInt aOff[],
                   const PylithInt aOff_x[],
                   const PylithScalar a[],
                   const PylithScalar a_t[],
                   const PylithScalar a_x[],
                   const PylithReal t,
                   const PylithReal s_tshift,
                   const PylithScalar x[],
                   const PylithReal n[],
                   const PylithInt numConstants,
                   const PylithScalar constants[],
                   PylithScalar Jf0[]) {
        assert(numS >= 1);
        assert(a);

        assert(numS >= 2);
        assert(Jf0);
        assert(sOff);

        const PylithInt spaceDim = dim+1; // :KLUDGE: dim passed in is spaceDim-1

        for (PylithInt i = 0; i < spaceDim; ++i) {
            Jf0[i*spaceDim+i] += +1.0;
        } // for
    } // Jf0ll_neg

    // ------------------------------------------------------------------------------------------------
    // Jf0 function for dynamic slip constraint equation for positive side of the fault.
    static inline
    void Jf0ll_pos(const PylithInt dim,
                   const PylithInt numS,
                   const PylithInt numA,
                   const PylithInt sOff[],
                   const PylithInt sOff_x[],
                   const PylithScalar s[],
                   const PylithScalar s_t[],
                   const PylithScalar s_x[],
                   const PylithInt aOff[],
                   const PylithInt aOff_x[],
                   const PylithScalar a[],
                   const PylithScalar a_t[],
                   const PylithScalar a_x[],
                   const PylithReal t,
                   const PylithReal s_tshift,
                   const PylithScalar x[],
                   const PylithReal n[],
                   const PylithInt numConstants,
                   const PylithScalar constants[],
                   PylithScalar Jf0[]) {
        assert(numS >= 1);
        assert(a);

        assert(numS >= 2);
        assert(Jf0);
        assert(sOff);

        const PylithInt spaceDim = dim; // dim passed is spatial dimension

        for (PylithInt i = 0; i < spaceDim; ++i) {
            Jf0[i*spaceDim+i] += +1.0;
        } // for
    } // Jf0ll_pos

    // Kernels for derived subfields

    // ------------------------------------------------------------------------------------------------
    /** Traction change from Lagrange multiplier field.
     *
     *  Solution fields: [disp(dim), vel(dim), lagrange(dim)]
     */
    static inline
    void tractionChange_asVector(const PylithInt dim,
                                 const PylithInt numS,
                                 const PylithInt numA,
                                 const PylithInt sOff[],
                                 const PylithInt sOff_x[],
                                 const PylithScalar s[],
                                 const PylithScalar s_t[],
                                 const PylithScalar s_x[],
                                 const PylithInt aOff[],
                                 const PylithInt aOff_x[],
                                 const PylithScalar a[],
                                 const PylithScalar a_t[],
                                 const PylithScalar a_x[],
                                 const PylithReal t,
                                 const PylithScalar x[],
                                 const PylithScalar n[],
                                 const PylithInt numConstants,
                                 const PylithScalar constants[],
                                 PylithScalar tractionChange[]) {
        const PylithInt spaceDim = dim; // dim passed in is embedded dimension (spaceDim).

        // Incoming solution fields.
        const PylithInt i_lagrange = numS-1;

        assert(numS >= 1);
        assert(sOff);
        assert(sOff[i_lagrange] >= 0);
        assert(tractionChange);

        const PylithScalar* lagrange = &s[sOff[i_lagrange]];

        switch (spaceDim) {
        case 2: {
            const PylithScalar tanDir[2] = {-n[1], n[0] };

            // 0=opening, 1=left-lateral
            tractionChange[0] = n[0]*lagrange[0] + n[1]*lagrange[1];
            tractionChange[1] = tanDir[0]*lagrange[0] + tanDir[1]*lagrange[1];
            break;
        } // case 2
        case 3: {
            const PylithScalar* refDir1 = &constants[0];
            const PylithScalar* refDir2 = &constants[3];
            PylithScalar tanDir1[3], tanDir2[3];
            pylith::fekernels::BoundaryDirections::tangential_directions(tanDir1, tanDir2, refDir1, refDir2, n);

            // 0=opening, 1=left-lateral, 2=reverse
            tractionChange[0] = n[0]*lagrange[0] + n[1]*lagrange[1] + n[2]*lagrange[2];
            tractionChange[1] = tanDir1[0]*lagrange[0] + tanDir1[1]*lagrange[1] + tanDir1[2]*lagrange[2];
            tractionChange[2] = tanDir2[0]*lagrange[0] + tanDir2[1]*lagrange[1] + tanDir2[2]*lagrange[2];
            break;
        } // case 3
        default:
            assert(0);
        } // switch

    } // tractionChange_asVector

}; // FaultCohesiveKin

/* End of file */
