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
 * Kernels for computing value from parameters for time dependent boundary conditions.
 *
 * f_0(x) +
 * \dot{f}_1(x) * H(t-t_1(x)) +
 * f_2(x) * a(t-t_2(x)) * H(t-t_2(s).
 */

#include "pylith/fekernels/fekernelsfwd.hh" // forward declarations

#include "pylith/utils/types.hh"

#include <cassert> // USES assert()

class pylith::fekernels::TimeDependentFn {
    // PUBLIC MEMBERS /////////////////////////////////////////////////////////////////////////////
public:

    // ============================================================================================
    // Boundary kernels
    // ============================================================================================

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
     * @param[in] n Face normal at point evaluation.
     * @param[in] numConstants Number of registered constants.
     * @param[in] constants Array of registered constants.
     * @param[out] value [dim] Output value.
     */

    // --------------------------------------------------------------------------------------------
    /** Scalar initial value term for time-dependent boundary condition.
     *
     * f_0(x)
     */
    static inline
    void initial_scalar_boundary(const PylithInt dim,
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
                                 PylithScalar value[]) {
        const PylithInt _numA = 1;
        assert(_numA == numA);
        assert(aOff);
        assert(a);

        const PylithInt i_initial = aOff[0];
        zero(1, value);
        initial_scalar_term(a[i_initial], value);
    } // initial_scalar

    // --------------------------------------------------------------------------------------------
    /** Vector initial value term for time-dependent boundary condition.
     *
     * f_0(x)
     */
    static inline
    void initial_vector_boundary(const PylithInt dim,
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
                                 PylithScalar value[]) {
        const PylithInt _numA = 1;
        assert(_numA == numA);
        assert(aOff);
        assert(a);

        const PylithInt i_initial = aOff[0];
        zero(dim, value);
        initial_vector_term(dim, &a[i_initial], value);
    } // initial_vector

    // --------------------------------------------------------------------------------------------
    /** Scalar rate term for time-dependent boundary condition.
     */
    static inline
    void rate_scalar_boundary(const PylithInt dim,
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
                              PylithScalar value[]) {
        const PylithInt _numA = 2;
        assert(_numA == numA);
        assert(aOff);
        assert(a);
        assert(value);

        const PylithInt i_rate = aOff[0];
        const PylithInt i_start = aOff[1];
        zero(1, value);
        rate_scalar_term(a[i_rate], a[i_start], t, value);
    } // rate_scalar

    // --------------------------------------------------------------------------------------------
    /** Vector rate term for time-dependent boundary condition.
     */
    static inline
    void rate_vector_boundary(const PylithInt dim,
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
                              PylithScalar value[]) {
        const PylithInt _numA = 2;
        assert(_numA == numA);
        assert(aOff);
        assert(a);

        const PylithInt i_rate = aOff[0];
        const PylithInt i_start = aOff[1];
        zero(dim, value);
        rate_vector_term(dim, &a[i_rate], a[i_start], t, value);
    } // rate_vector

    // --------------------------------------------------------------------------------------------
    /** Scalar time history term for time-dependent boundary condition.
     */
    static inline
    void timeHistory_scalar_boundary(const PylithInt dim,
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
                                     PylithScalar value[]) {
        const PylithInt _numA = 3;
        assert(_numA == numA);
        assert(aOff);
        assert(a);

        const PylithInt i_amplitude = aOff[0];
        const PylithInt i_value = aOff[2];
        zero(1, value);
        timeHistory_scalar_term(a[i_amplitude], a[i_value], value);
    } // timeHistory_scalar

    // --------------------------------------------------------------------------------------------
    /** Vector time history term for time-dependent boundary condition.
     */
    static inline
    void timeHistory_vector_boundary(const PylithInt dim,
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
                                     PylithScalar value[]) {
        const PylithInt _numA = 3;
        assert(_numA == numA);
        assert(aOff);
        assert(a);

        const PylithInt i_amplitude = aOff[0];
        const PylithInt i_value = aOff[2];
        zero(dim, value);
        timeHistory_vector_term(dim, &a[i_amplitude], a[i_value], value);
    } // timeHistory_vector

    // --------------------------------------------------------------------------------------------
    /** Compute boundary condition scalar value using initial and rate terms.
     *
     * value = initial_term + rate_term
     */
    static inline
    void initialRate_scalar_boundary(const PylithInt dim,
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
                                     PylithScalar value[]) {
        const PylithInt _numA = 3;
        assert(_numA == numA);
        assert(aOff);
        assert(a);
        assert(value);

        const PylithInt i_initial = aOff[0];
        const PylithInt i_rate = aOff[1];
        const PylithInt i_start = aOff[2];
        zero(1, value);
        initial_scalar_term(a[i_initial], value);
        rate_scalar_term(a[i_rate], a[i_start], t, value);
    } // initialRate_scalar

    // --------------------------------------------------------------------------------------------
    /** Compute boundary condition vector value using initial and rate terms.
     *
     * value = initial_term + rate_term
     */
    static inline
    void initialRate_vector_boundary(const PylithInt dim,
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
                                     PylithScalar value[]) {
        const PylithInt _numA = 3;
        assert(_numA == numA);
        assert(aOff);
        assert(a);

        const PylithInt i_initial = aOff[0];
        const PylithInt i_rate = aOff[1];
        const PylithInt i_start = aOff[2];
        zero(dim, value);
        initial_vector_term(dim, &a[i_initial], value);
        rate_vector_term(dim, &a[i_rate], a[i_start], t, value);
    } // initialRate_vector

    // --------------------------------------------------------------------------------------------
    /** Compute boundary condition scalar value using initial and time history terms.
     *
     * value = inital_term + timeHistory+term
     */
    static inline
    void initialTimeHistory_scalar_boundary(const PylithInt dim,
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
                                            PylithScalar value[]) {
        const PylithInt _numA = 4;
        assert(_numA == numA);
        assert(aOff);
        assert(a);

        const PylithInt i_initial = aOff[0];
        const PylithInt i_thAmp = aOff[1];
        const PylithInt i_thValue = aOff[3];
        zero(1, value);
        initial_scalar_term(a[i_initial], value);
        timeHistory_scalar_term(a[i_thAmp], a[i_thValue], value);
    } // initialTimeHistory_scalar

    // --------------------------------------------------------------------------------------------
    /** Compute boundary condition vector value using initial and time history terms.
     *
     * f_0(x) +
     * \dot{f}_1(x) * H(t-t_1(x)) +
     * f_2(x) * a(t-t_2(x)) * H(t-t_2(s).
     */
    static inline
    void initialTimeHistory_vector_boundary(const PylithInt dim,
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
                                            PylithScalar value[]) {
        const PylithInt _numA = 4;
        assert(_numA == numA);
        assert(aOff);
        assert(a);
        assert(value);

        const PylithInt i_initial = aOff[0];
        const PylithInt i_thAmp = aOff[1];
        const PylithInt i_thValue = aOff[3];
        zero(dim, value);
        initial_vector_term(dim, &a[i_initial], value);
        timeHistory_vector_term(dim, &a[i_thAmp], a[i_thValue], value);
    } // initialTimeHistory_vector

    // --------------------------------------------------------------------------------------------
    /** Compute boundary condition scalar value using rate and time history terms.
     *
     * value = rate_term + timeHistory_term
     */
    static inline
    void rateTimeHistory_scalar_boundary(const PylithInt dim,
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
                                         PylithScalar value[]) {
        const PylithInt _numA = 5;
        assert(_numA == numA);
        assert(aOff);
        assert(a);

        const PylithInt i_rateAmp = aOff[0];
        const PylithInt i_rateStart = aOff[1];
        const PylithInt i_thAmp = aOff[2];
        const PylithInt i_thValue = aOff[4];
        zero(1, value);
        rate_scalar_term(a[i_rateAmp], a[i_rateStart], t, value);
        timeHistory_scalar_term(a[i_thAmp], a[i_thValue], value);
    } // rateTimeHistory_scalar

    // --------------------------------------------------------------------------------------------
    /** Compute boundary condition vector value using rate and time history terms.
     *
     * value = rate_term + timeHistory_term
     */
    static inline
    void rateTimeHistory_vector_boundary(const PylithInt dim,
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
                                         PylithScalar value[]) {
        const PylithInt _numA = 5;
        assert(_numA == numA);
        assert(aOff);
        assert(a);

        const PylithInt i_rateAmp = aOff[0];
        const PylithInt i_rateStart = aOff[1];
        const PylithInt i_thAmp = aOff[2];
        const PylithInt i_thValue = aOff[4];
        zero(dim, value);
        rate_vector_term(dim, &a[i_rateAmp], a[i_rateStart], t, value);
        timeHistory_vector_term(dim, &a[i_thAmp], a[i_thValue], value);
    } // rateTimeHistory_vector

    // --------------------------------------------------------------------------------------------
    /** Compute boundary condition scalar value using initial, rate ,and time history terms.
     *
     * value = initial_term + rate_term + timeHistory_term
     */
    static inline
    void initialRateTimeHistory_scalar_boundary(const PylithInt dim,
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
                                                PylithScalar value[]) {
        const PylithInt _numA = 6;
        assert(_numA == numA);
        assert(aOff);
        assert(a);

        const PylithInt i_initialAmp = aOff[0];
        const PylithInt i_rateAmp = aOff[1];
        const PylithInt i_rateStart = aOff[2];
        const PylithInt i_thAmp = aOff[3];
        const PylithInt i_thValue = aOff[5];
        zero(1, value);
        initial_scalar_term(a[i_initialAmp], value);
        rate_scalar_term(a[i_rateAmp], a[i_rateStart], t, value);
        timeHistory_scalar_term(a[i_thAmp], a[i_thValue], value);
    } // initialRateTimeHistory_scalar

    // --------------------------------------------------------------------------------------------
    /** Compute boundary condition vector value using initial, rate, and time history terms.
     *
     * value = initial_term + rate_term + timeHistory_term
     */
    static inline
    void initialRateTimeHistory_vector_boundary(const PylithInt dim,
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
                                                PylithScalar value[]) {
        const PylithInt _numA = 6;
        assert(_numA == numA);
        assert(aOff);
        assert(a);

        const PylithInt i_initialAmp = aOff[0];
        const PylithInt i_rateAmp = aOff[1];
        const PylithInt i_rateStart = aOff[2];
        const PylithInt i_thAmp = aOff[3];
        const PylithInt i_thValue = aOff[5];
        zero(dim, value);
        initial_vector_term(dim, &a[i_initialAmp], value);
        rate_vector_term(dim, &a[i_rateAmp], a[i_rateStart], t, value);
        timeHistory_vector_term(dim, &a[i_thAmp], a[i_thValue], value);
    } // initialRateTimeHistory_vector

    // ============================================================================================
    // Non-boundary kernels
    // ============================================================================================

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
     * @param[in] numConstants Number of registered constants.
     * @param[in] constants Array of registered constants.
     * @param[out] value [dim] Output value.
     */

    // --------------------------------------------------------------------------------------------
    /** Scalar initial value term for time-dependent boundary condition.
     *
     * f_0(x)
     */
    static inline
    void initial_scalar(const PylithInt dim,
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
                        const PylithInt numConstants,
                        const PylithScalar constants[],
                        PylithScalar value[]) {
        const PylithInt _numA = 1;
        assert(_numA == numA);
        assert(aOff);
        assert(a);

        const PylithInt i_initial = aOff[0];
        zero(1, value);
        initial_scalar_term(a[i_initial], value);
    } // initial_scalar

    // --------------------------------------------------------------------------------------------
    /** Vector initial value term for time-dependent boundary condition.
     *
     * f_0(x)
     */
    static inline
    void initial_vector(const PylithInt dim,
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
                        const PylithInt numConstants,
                        const PylithScalar constants[],
                        PylithScalar value[]) {
        const PylithInt _numA = 1;
        assert(_numA == numA);
        assert(aOff);
        assert(a);

        const PylithInt i_initial = aOff[0];
        zero(dim, value);
        initial_vector_term(dim, &a[i_initial], value);
    } // initial_vector

    // --------------------------------------------------------------------------------------------
    /** Scalar rate term for time-dependent boundary condition.
     */
    static inline
    void rate_scalar(const PylithInt dim,
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
                     const PylithInt numConstants,
                     const PylithScalar constants[],
                     PylithScalar value[]) {
        const PylithInt _numA = 2;
        assert(_numA == numA);
        assert(aOff);
        assert(a);
        assert(value);

        const PylithInt i_rate = aOff[0];
        const PylithInt i_start = aOff[1];
        zero(1, value);
        rate_scalar_term(a[i_rate], a[i_start], t, value);
    } // rate_scalar

    // --------------------------------------------------------------------------------------------
    /** Vector rate term for time-dependent boundary condition.
     */
    static inline
    void rate_vector(const PylithInt dim,
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
                     const PylithInt numConstants,
                     const PylithScalar constants[],
                     PylithScalar value[]) {
        const PylithInt _numA = 2;
        assert(_numA == numA);
        assert(aOff);
        assert(a);

        const PylithInt i_rate = aOff[0];
        const PylithInt i_start = aOff[1];
        zero(dim, value);
        rate_vector_term(dim, &a[i_rate], a[i_start], t, value);
    } // rate_vector

    // --------------------------------------------------------------------------------------------
    /** Scalar time history term for time-dependent boundary condition.
     */
    static inline
    void timeHistory_scalar(const PylithInt dim,
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
                            const PylithInt numConstants,
                            const PylithScalar constants[],
                            PylithScalar value[]) {
        const PylithInt _numA = 3;
        assert(_numA == numA);
        assert(aOff);
        assert(a);

        const PylithInt i_amplitude = aOff[0];
        const PylithInt i_value = aOff[2];
        zero(1, value);
        timeHistory_scalar_term(a[i_amplitude], a[i_value], value);
    } // timeHistory_scalar

    // --------------------------------------------------------------------------------------------
    /** Vector time history term for time-dependent boundary condition.
     */
    static inline
    void timeHistory_vector(const PylithInt dim,
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
                            const PylithInt numConstants,
                            const PylithScalar constants[],
                            PylithScalar value[]) {
        const PylithInt _numA = 3;
        assert(_numA == numA);
        assert(aOff);
        assert(a);

        const PylithInt i_amplitude = aOff[0];
        const PylithInt i_value = aOff[2];
        zero(dim, value);
        timeHistory_vector_term(dim, &a[i_amplitude], a[i_value], value);
    } // timeHistory_vector

    // --------------------------------------------------------------------------------------------
    /** Compute boundary condition scalar value using initial and rate terms.
     *
     * value = initial_term + rate_term
     */
    static inline
    void initialRate_scalar(const PylithInt dim,
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
                            const PylithInt numConstants,
                            const PylithScalar constants[],
                            PylithScalar value[]) {
        const PylithInt _numA = 3;
        assert(_numA == numA);
        assert(aOff);
        assert(a);
        assert(value);

        const PylithInt i_initial = aOff[0];
        const PylithInt i_rate = aOff[1];
        const PylithInt i_start = aOff[2];
        zero(1, value);
        initial_scalar_term(a[i_initial], value);
        rate_scalar_term(a[i_rate], a[i_start], t, value);
    } // initialRate_scalar

    // --------------------------------------------------------------------------------------------
    /** Compute boundary condition vector value using initial and rate terms.
     *
     * value = initial_term + rate_term
     */
    static inline
    void initialRate_vector(const PylithInt dim,
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
                            const PylithInt numConstants,
                            const PylithScalar constants[],
                            PylithScalar value[]) {
        const PylithInt _numA = 3;
        assert(_numA == numA);
        assert(aOff);
        assert(a);

        const PylithInt i_initial = aOff[0];
        const PylithInt i_rate = aOff[1];
        const PylithInt i_start = aOff[2];
        zero(dim, value);
        initial_vector_term(dim, &a[i_initial], value);
        rate_vector_term(dim, &a[i_rate], a[i_start], t, value);
    } // initialRate_vector

    // --------------------------------------------------------------------------------------------
    /** Compute boundary condition scalar value using initial and time history terms.
     *
     * value = inital_term + timeHistory+term
     */
    static inline
    void initialTimeHistory_scalar(const PylithInt dim,
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
                                   const PylithInt numConstants,
                                   const PylithScalar constants[],
                                   PylithScalar value[]) {
        const PylithInt _numA = 4;
        assert(_numA == numA);
        assert(aOff);
        assert(a);

        const PylithInt i_initial = aOff[0];
        const PylithInt i_thAmp = aOff[1];
        const PylithInt i_thValue = aOff[3];
        zero(1, value);
        initial_scalar_term(a[i_initial], value);
        timeHistory_scalar_term(a[i_thAmp], a[i_thValue], value);
    } // initialTimeHistory_scalar

    // --------------------------------------------------------------------------------------------
    /** Compute boundary condition vector value using initial and time history terms.
     *
     * f_0(x) +
     * \dot{f}_1(x) * H(t-t_1(x)) +
     * f_2(x) * a(t-t_2(x)) * H(t-t_2(s).
     */
    static inline
    void initialTimeHistory_vector(const PylithInt dim,
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
                                   const PylithInt numConstants,
                                   const PylithScalar constants[],
                                   PylithScalar value[]) {
        const PylithInt _numA = 4;
        assert(_numA == numA);
        assert(aOff);
        assert(a);
        assert(value);

        const PylithInt i_initial = aOff[0];
        const PylithInt i_thAmp = aOff[1];
        const PylithInt i_thValue = aOff[3];
        zero(dim, value);
        initial_vector_term(dim, &a[i_initial], value);
        timeHistory_vector_term(dim, &a[i_thAmp], a[i_thValue], value);
    } // initialTimeHistory_vector

    // --------------------------------------------------------------------------------------------
    /** Compute boundary condition scalar value using rate and time history terms.
     *
     * value = rate_term + timeHistory_term
     */
    static inline
    void rateTimeHistory_scalar(const PylithInt dim,
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
                                const PylithInt numConstants,
                                const PylithScalar constants[],
                                PylithScalar value[]) {
        const PylithInt _numA = 5;
        assert(_numA == numA);
        assert(aOff);
        assert(a);

        const PylithInt i_rateAmp = aOff[0];
        const PylithInt i_rateStart = aOff[1];
        const PylithInt i_thAmp = aOff[2];
        const PylithInt i_thValue = aOff[4];
        zero(1, value);
        rate_scalar_term(a[i_rateAmp], a[i_rateStart], t, value);
        timeHistory_scalar_term(a[i_thAmp], a[i_thValue], value);
    } // rateTimeHistory_scalar

    // --------------------------------------------------------------------------------------------
    /** Compute boundary condition vector value using rate and time history terms.
     *
     * value = rate_term + timeHistory_term
     */
    static inline
    void rateTimeHistory_vector(const PylithInt dim,
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
                                const PylithInt numConstants,
                                const PylithScalar constants[],
                                PylithScalar value[]) {
        const PylithInt _numA = 5;
        assert(_numA == numA);
        assert(aOff);
        assert(a);

        const PylithInt i_rateAmp = aOff[0];
        const PylithInt i_rateStart = aOff[1];
        const PylithInt i_thAmp = aOff[2];
        const PylithInt i_thValue = aOff[4];
        zero(dim, value);
        rate_vector_term(dim, &a[i_rateAmp], a[i_rateStart], t, value);
        timeHistory_vector_term(dim, &a[i_thAmp], a[i_thValue], value);
    } // rateTimeHistory_vector

    // --------------------------------------------------------------------------------------------
    /** Compute boundary condition scalar value using initial, rate ,and time history terms.
     *
     * value = initial_term + rate_term + timeHistory_term
     */
    static inline
    void initialRateTimeHistory_scalar(const PylithInt dim,
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
                                       const PylithInt numConstants,
                                       const PylithScalar constants[],
                                       PylithScalar value[]) {
        const PylithInt _numA = 6;
        assert(_numA == numA);
        assert(aOff);
        assert(a);

        const PylithInt i_initialAmp = aOff[0];
        const PylithInt i_rateAmp = aOff[1];
        const PylithInt i_rateStart = aOff[2];
        const PylithInt i_thAmp = aOff[3];
        const PylithInt i_thValue = aOff[5];
        zero(1, value);
        initial_scalar_term(a[i_initialAmp], value);
        rate_scalar_term(a[i_rateAmp], a[i_rateStart], t, value);
        timeHistory_scalar_term(a[i_thAmp], a[i_thValue], value);
    } // initialRateTimeHistory_scalar

    // --------------------------------------------------------------------------------------------
    /** Compute boundary condition vector value using initial, rate, and time history terms.
     *
     * value = initial_term + rate_term + timeHistory_term
     */
    static inline
    void initialRateTimeHistory_vector(const PylithInt dim,
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
                                       const PylithInt numConstants,
                                       const PylithScalar constants[],
                                       PylithScalar value[]) {
        const PylithInt _numA = 6;
        assert(_numA == numA);
        assert(aOff);
        assert(a);

        const PylithInt i_initialAmp = aOff[0];
        const PylithInt i_rateAmp = aOff[1];
        const PylithInt i_rateStart = aOff[2];
        const PylithInt i_thAmp = aOff[3];
        const PylithInt i_thValue = aOff[5];
        zero(dim, value);
        initial_vector_term(dim, &a[i_initialAmp], value);
        rate_vector_term(dim, &a[i_rateAmp], a[i_rateStart], t, value);
        timeHistory_vector_term(dim, &a[i_thAmp], a[i_thValue], value);
    } // initialRateTimeHistory_vector

    // ============================================================================================
    // Time-dependent terms
    // ============================================================================================

    // --------------------------------------------------------------------------------------------
    /** Scalar initial value term.
     *
     * value = f_0(x)
     */
    static inline
    void initial_scalar_term(const PylithReal initialValue,
                             PylithScalar value[]) {
        assert(value);

        value[0] += initialValue;
    } // initial_scalar_term

    // --------------------------------------------------------------------------------------------
    /** Vector initial value term.
     *
     * value = f_0(x)
     */
    static inline
    void initial_vector_term(const PylithInt dim,
                             const PylithReal initialValue[],
                             PylithScalar value[]) {
        assert(initialValue);
        assert(value);

        for (PylithInt i = 0; i < dim; ++i) {
            value[i] += initialValue[i];
        } // for
    } // initial_vector_term

    // --------------------------------------------------------------------------------------------
    /** Scalar rate term for time-dependent boundary condition.
     *
     * value = \dot{f}_1(x) * (t-t_1(x)) for t >= t_1(x).
     */
    static inline
    void rate_scalar_term(const PylithReal rateValue,
                          const PylithReal startTime,
                          const PylithReal t,
                          PylithScalar value[]) {
        assert(value);

        const PylithScalar tRel = t - startTime;
        if (tRel > 0.0) {
            value[0] += rateValue*tRel;
        } // if
    } // rate_scalar_term

    // --------------------------------------------------------------------------------------------
    /** Vector rate term for time-dependent boundary condition.
     *
     * value = \dot{f}_1(x) * (t-t_1(x)) for t >= t_1(x).
     */
    static inline
    void rate_vector_term(const PylithInt dim,
                          const PylithReal rateValue[],
                          const PylithReal startTime,
                          const PylithReal t,
                          PylithScalar value[]) {
        assert(rateValue);
        assert(value);

        const PylithScalar tRel = t - startTime;
        if (tRel > 0.0) {
            for (PylithInt i = 0; i < dim; ++i) {
                value[i] += rateValue[i] * tRel;
            } // for
        } // if
    } // rate_vector_term

    // --------------------------------------------------------------------------------------------
    /** Scalar time history term for time-dependent boundary condition.
     *
     * value = f_2(x) * a(t-t_2(x)) for t >= t_2(x).
     */
    static inline
    void timeHistory_scalar_term(const PylithReal amplitude,
                                 const PylithReal timeHistoryValue,
                                 PylithScalar value[]) {
        assert(value);

        value[0] += amplitude * timeHistoryValue;
    } // timeHistory_scalar_term

    // --------------------------------------------------------------------------------------------
    /** Vector time history term for time-dependent boundary condition.
     *
     * value = f_2(x) * a(t-t_2(x)) for t >= t_2(x).
     */
    static inline
    void timeHistory_vector_term(const PylithInt dim,
                                 const PylithReal amplitude[],
                                 const PylithReal timeHistoryValue,
                                 PylithScalar value[]) {
        assert(amplitude);
        assert(value);

        for (PylithInt i = 0; i < dim; ++i) {
            value[i] += amplitude[i] * timeHistoryValue;
        } // for
    } // timeHistory_vector_term

    // --------------------------------------------------------------------------------------------
    // Set scalar value to zero.
    static inline
    void zero(const PylithInt dim,
              PylithScalar value[]) {
        assert(value);

        for (PylithInt i = 0; i < dim; ++i) {
            value[i] = 0.0;
        } // for
    } // zero

}; // TimeDependentFn

// End of file
