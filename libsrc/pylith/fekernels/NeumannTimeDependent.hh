/*
 * ================================================================================================
 * This code is part of PyLith, developed through the Computational Infrastructure
 * for Geodynamics (https://github.com/geodynamics/pylith).
 *
 * Copyright (c) 2010-2024, University of California, Davis and the PyLith Development Team.
 * All rights reserved.
 *
 * See https://mit-license.org/ and LICENSE.md and for license information.
 * =================================================================================================
 */
#pragma once

/** @file libsrc/fekernels/NeumannTimeDependent.h
 *
 * Kernels for computing value from parameters for Neumann time dependent boundary conditions.
 *
 * \int_{\Gamma_\tau} \trialvec[u] \vec{\tau}(\vec{x},t) d\Gamma
 */

// Include directives ---------------------------------------------------
#include "pylith/fekernels/fekernelsfwd.hh" // forward declarations

#include "pylith/fekernels/BoundaryDirections.hh" // USES tangential_directions()
#include "pylith/fekernels/TimeDependentFn.hh" // USES TimeDependentFn

#include "pylith/utils/types.hh"

#include <cassert> // USES assert()

class pylith::fekernels::NeumannTimeDependent {
public:

    /** Kernel interface.
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
     * @param[in] n Unit vector normal to boundary.
     * @param[in] numConstants Number of registered constants.
     * @param[in] constants Array of registered constants.
     * @param[out] f0 [dim].
     */

    // --------------------------------------------------------------------------------------------
    /** Scalar initial value term for time-dependent boundary condition.
     *
     * f_0(x)
     */
    static inline
    void f0_initial_scalar(const PylithInt dim,
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
                           const PylithReal x[],
                           const PylithReal n[],
                           const PylithInt numConstants,
                           const PylithScalar constants[],
                           PylithScalar f0[]) {
        pylith::fekernels::TimeDependentFn::initial_scalar_boundary(
            dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, n, numConstants, constants, f0);
    } // f0_initial_scalar

    // --------------------------------------------------------------------------------------------
    /** Vector initial value term for time-dependent boundary condition.
     *
     * f_0(x)
     */
    static inline
    void f0_initial_vector(const PylithInt dim,
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
                           const PylithReal x[],
                           const PylithReal n[],
                           const PylithInt numConstants,
                           const PylithScalar constants[],
                           PylithScalar f0[]) {
        PylithScalar values[3] = { 0.0, 0.0, 0.0 };
        pylith::fekernels::TimeDependentFn::initial_vector_boundary(
            dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, n, numConstants, constants, values);

        switch (dim) {
        case 2:
            pylith::fekernels::BoundaryDirections::toXY(f0, values, n);
            break;
        case 3: {
            const PylithScalar* refDir1 = &constants[0];
            const PylithScalar* refDir2 = &constants[3];
            pylith::fekernels::BoundaryDirections::toXYZ(f0, values, refDir1, refDir2, n);
            break;
        }
        default:
            assert(0);
        } // switch
    } // f0_initial_vector

    // --------------------------------------------------------------------------------------------
    /** Scalar rate term for time-dependent boundary condition.
     *
     * \dot{f}_1(x) * (t-t_1(x)) for t >= t_1(x).
     */
    static inline
    void f0_rate_scalar(const PylithInt dim,
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
                        const PylithReal x[],
                        const PylithReal n[],
                        const PylithInt numConstants,
                        const PylithScalar constants[],
                        PylithScalar f0[]) {
        pylith::fekernels::TimeDependentFn::rate_scalar_boundary(
            dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, n, numConstants, constants, f0);
    } // f0_rate_scalar

    // --------------------------------------------------------------------------------------------
    /** Vector rate term for time-dependent boundary condition.
     *
     * \dot{f}_1(x) * (t-t_1(x)) for t >= t_1(x).
     */
    static inline
    void f0_rate_vector(const PylithInt dim,
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
                        const PylithReal x[],
                        const PylithReal n[],
                        const PylithInt numConstants,
                        const PylithScalar constants[],
                        PylithScalar f0[]) {
        PylithScalar values[3] = { 0.0, 0.0, 0.0 };
        pylith::fekernels::TimeDependentFn::rate_vector_boundary(
            dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, n, numConstants, constants, values);

        switch (dim) {
        case 2:
            pylith::fekernels::BoundaryDirections::toXY(f0, values, n);
            break;
        case 3: {
            const PylithScalar* refDir1 = &constants[0];
            const PylithScalar* refDir2 = &constants[3];
            pylith::fekernels::BoundaryDirections::toXYZ(f0, values, refDir1, refDir2, n);
            break;
        }
        default:
            assert(0);
        } // switch
    } // f0_rate_vector

    // --------------------------------------------------------------------------------------------
    /** Scalar time history term for time-dependent boundary condition.
     *
     * f_2(x) * a(t-t_2(x)) for t >= t_2(x).
     */
    static inline
    void f0_timeHistory_scalar(const PylithInt dim,
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
                               const PylithReal x[],
                               const PylithReal n[],
                               const PylithInt numConstants,
                               const PylithScalar constants[],
                               PylithScalar f0[]) {
        pylith::fekernels::TimeDependentFn::timeHistory_scalar_boundary(
            dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, n, numConstants, constants, f0);
    } // f0_timeHistory_scalar

    // --------------------------------------------------------------------------------------------
    /** Vector time history term for time-dependent boundary condition.
     *
     * f_2(x) * a(t-t_2(x)) for t >= t_2(x).
     */
    static inline
    void f0_timeHistory_vector(const PylithInt dim,
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
                               const PylithReal x[],
                               const PylithReal n[],
                               const PylithInt numConstants,
                               const PylithScalar constants[],
                               PylithScalar f0[]) {
        PylithScalar values[3] = { 0.0, 0.0, 0.0 };
        pylith::fekernels::TimeDependentFn::timeHistory_vector_boundary(
            dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, n, numConstants, constants, values);

        switch (dim) {
        case 2:
            pylith::fekernels::BoundaryDirections::toXY(f0, values, n);
            break;
        case 3: {
            const PylithScalar* refDir1 = &constants[0];
            const PylithScalar* refDir2 = &constants[3];
            pylith::fekernels::BoundaryDirections::toXYZ(f0, values, refDir1, refDir2, n);
            break;
        }
        default:
            assert(0);
        } // switch
    } // f0_timeHistory_vector

    // --------------------------------------------------------------------------------------------
    /** Compute boundary condition scalar value using initial and rate terms.
     *
     * f_0(x) +
     * \dot{f}_1(x) * H(t-t_1(x)) +
     */
    static inline
    void f0_initialRate_scalar(const PylithInt dim,
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
                               const PylithReal x[],
                               const PylithReal n[],
                               const PylithInt numConstants,
                               const PylithScalar constants[],
                               PylithScalar f0[]) {
        pylith::fekernels::TimeDependentFn::initialRate_scalar_boundary(
            dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, n, numConstants, constants, f0);
    } // f0_initialRate_scalar

    // --------------------------------------------------------------------------------------------
    /** Compute boundary condition vector value using initial and rate terms.
     *
     * f_0(x) +
     * \dot{f}_1(x) * H(t-t_1(x)) +
     */
    static inline
    void f0_initialRate_vector(const PylithInt dim,
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
                               const PylithReal x[],
                               const PylithReal n[],
                               const PylithInt numConstants,
                               const PylithScalar constants[],
                               PylithScalar f0[]) {
        PylithScalar values[3] = { 0.0, 0.0, 0.0 };
        pylith::fekernels::TimeDependentFn::initialRate_vector_boundary(
            dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, n, numConstants, constants, values);

        switch (dim) {
        case 2:
            pylith::fekernels::BoundaryDirections::toXY(f0, values, n);
            break;
        case 3: {
            const PylithScalar* refDir1 = &constants[0];
            const PylithScalar* refDir2 = &constants[3];
            pylith::fekernels::BoundaryDirections::toXYZ(f0, values, refDir1, refDir2, n);
            break;
        }
        default:
            assert(0);
        } // switch
    } // f0_initialRate_vector

    // --------------------------------------------------------------------------------------------
    /** Compute boundary condition scalar value using initial and time history terms.
     *
     * f_0(x) +
     * \dot{f}_1(x) * H(t-t_1(x)) +
     * f_2(x) * a(t-t_2(x)) * H(t-t_2(s).
     */
    static inline
    void f0_initialTimeHistory_scalar(const PylithInt dim,
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
                                      const PylithReal x[],
                                      const PylithReal n[],
                                      const PylithInt numConstants,
                                      const PylithScalar constants[],
                                      PylithScalar f0[]) {
        pylith::fekernels::TimeDependentFn::initialTimeHistory_scalar_boundary(
            dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, n, numConstants, constants, f0);
    } // f0_initialTimeHistory_scalar

    // --------------------------------------------------------------------------------------------
    /** Compute boundary condition vector value using initial and time history terms.
     *
     * f_0(x) +
     * \dot{f}_1(x) * H(t-t_1(x)) +
     * f_2(x) * a(t-t_2(x)) * H(t-t_2(s).
     */
    static inline
    void f0_initialTimeHistory_vector(const PylithInt dim,
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
                                      const PylithReal x[],
                                      const PylithReal n[],
                                      const PylithInt numConstants,
                                      const PylithScalar constants[],
                                      PylithScalar f0[]) {
        PylithScalar values[3] = { 0.0, 0.0, 0.0 };
        pylith::fekernels::TimeDependentFn::initialTimeHistory_vector_boundary(
            dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, n, numConstants, constants, values);

        switch (dim) {
        case 2:
            pylith::fekernels::BoundaryDirections::toXY(f0, values, n);
            break;
        case 3: {
            const PylithScalar* refDir1 = &constants[0];
            const PylithScalar* refDir2 = &constants[3];
            pylith::fekernels::BoundaryDirections::toXYZ(f0, values, refDir1, refDir2, n);
            break;
        }
        default:
            assert(0);
        } // switch
    } // f0_initialTimeHistory_vector

    // --------------------------------------------------------------------------------------------
    /** Compute boundary condition scalar value using rate and time history terms.
     *
     * \dot{f}_1(x) * H(t-t_1(x)) +
     * f_2(x) * a(t-t_2(x)) * H(t-t_2(s).
     */
    static inline
    void f0_rateTimeHistory_scalar(const PylithInt dim,
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
                                   const PylithReal x[],
                                   const PylithReal n[],
                                   const PylithInt numConstants,
                                   const PylithScalar constants[],
                                   PylithScalar f0[]) {
        pylith::fekernels::TimeDependentFn::rateTimeHistory_scalar_boundary(
            dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, n, numConstants, constants, f0);
    } // f0_rateTimeHistory_scalar

    // --------------------------------------------------------------------------------------------
    /** Compute boundary condition vector value using rate and time history terms.
     *
     * \dot{f}_1(x) * H(t-t_1(x)) +
     * f_2(x) * a(t-t_2(x)) * H(t-t_2(s).
     */
    static inline
    void f0_rateTimeHistory_vector(const PylithInt dim,
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
                                   const PylithReal x[],
                                   const PylithReal n[],
                                   const PylithInt numConstants,
                                   const PylithScalar constants[],
                                   PylithScalar f0[]) {
        PylithScalar values[3] = { 0.0, 0.0, 0.0 };
        pylith::fekernels::TimeDependentFn::rateTimeHistory_vector_boundary(
            dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, n, numConstants, constants, values);

        switch (dim) {
        case 2:
            pylith::fekernels::BoundaryDirections::toXY(f0, values, n);
            break;
        case 3: {
            const PylithScalar* refDir1 = &constants[0];
            const PylithScalar* refDir2 = &constants[3];
            pylith::fekernels::BoundaryDirections::toXYZ(f0, values, refDir1, refDir2, n);
            break;
        }
        default:
            assert(0);
        } // switch
    } // f0_rateTimeHistory_vector

    // --------------------------------------------------------------------------------------------
    /** Compute boundary condition scalar value using initial, rate ,and time history terms.
     *
     * f_0(x) +
     * \dot{f}_1(x) * H(t-t_1(x)) +
     * f_2(x) * a(t-t_2(x)) * H(t-t_2(s).
     */
    static inline
    void f0_initialRateTimeHistory_scalar(const PylithInt dim,
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
                                          const PylithReal x[],
                                          const PylithReal n[],
                                          const PylithInt numConstants,
                                          const PylithScalar constants[],
                                          PylithScalar f0[]) {
        pylith::fekernels::TimeDependentFn::initialRateTimeHistory_scalar_boundary(
            dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, n, numConstants, constants, f0);
    } // f0_initialRateTimeHistory_scalar

    // --------------------------------------------------------------------------------------------
    /** Compute boundary condition vector value using initial, rate, and time history terms.
     *
     * f_0(x) +
     * \dot{f}_1(x) * H(t-t_1(x)) +
     * f_2(x) * a(t-t_2(x)) * H(t-t_2(s).
     */
    static inline
    void f0_initialRateTimeHistory_vector(const PylithInt dim,
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
                                          const PylithReal x[],
                                          const PylithReal n[],
                                          const PylithInt numConstants,
                                          const PylithScalar constants[],
                                          PylithScalar f0[]) {
        PylithScalar values[3] = { 0.0, 0.0, 0.0 };
        pylith::fekernels::TimeDependentFn::initialRateTimeHistory_vector_boundary(
            dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, n, numConstants, constants, values);

        switch (dim) {
        case 2:
            pylith::fekernels::BoundaryDirections::toXY(f0, values, n);
            break;
        case 3: {
            const PylithScalar* refDir1 = &constants[0];
            const PylithScalar* refDir2 = &constants[3];
            pylith::fekernels::BoundaryDirections::toXYZ(f0, values, refDir1, refDir2, n);
            break;
        }
        default:
            assert(0);
        } // switch
    } // f0_initialRateTimeHistory_vector

}; // NeumannTimeDependent

// End of file
