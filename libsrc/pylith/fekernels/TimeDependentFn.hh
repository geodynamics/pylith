/* -*- C -*-
 *
 * ----------------------------------------------------------------------
 *
 * Brad T. Aagaard, U.S. Geological Survey
 * Charles A. Williams, GNS Science
 * Matthew G. Knepley, University at Buffalo
 *
 * This code was developed as part of the Computational Infrastructure
 * for Geodynamics (http:*geodynamics.org).
 *
 * Copyright (c) 2010-2022 University of California, Davis
 *
 * See LICENSE.md for license information.
 *
 * ----------------------------------------------------------------------
 */

/** @file libsrc/fekernels/TimeDependentFn.hh
 *
 * Kernels for computing value from parameters for time dependent boundary conditions.
 *
 * f_0(x) +
 * \dot{f}_1(x) * H(t-t_1(x)) +
 * f_2(x) * a(t-t_2(x)) * H(t-t_2(s).
 */

#if !defined(pylith_fekernels_timedependentfn_hh)
#define pylith_fekernels_timedependentfn_hh

// Include directives ---------------------------------------------------
#include "fekernelsfwd.hh" // forward declarations

#include "pylith/utils/types.hh"

#include <cassert> // USES assert()

class pylith::fekernels::TimeDependentFn {
    // PUBLIC MEMBERS ///////////////////////////////////////////////////////
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
     * @param[in] n Face normal at point evaluation.
     * @param[in] numConstants Number of registered constants.
     * @param[in] constants Array of registered constants.
     * @param[out] f0 [dim].
     */

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
                        const PylithScalar n[],
                        const PylithInt numConstants,
                        const PylithScalar constants[],
                        PylithScalar value[]) {
        const PylithInt _numA = 1;
        assert(_numA == numA);
        assert(aOff);
        assert(a);
        assert(value);

        const PylithInt i_initial = aOff[0];

        value[0] = a[i_initial];
    } // initial_scalar

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
                        const PylithScalar n[],
                        const PylithInt numConstants,
                        const PylithScalar constants[],
                        PylithScalar value[]) {
        const PylithInt _numA = 1;
        assert(_numA == numA);
        assert(aOff);
        assert(a);
        assert(value);

        const PylithInt i_initial = aOff[0];

        for (PylithInt i = 0; i < dim; ++i) {
            value[i] = a[i_initial+i];
        } // for
    } // initial_vector

    /** Scalar rate term for time-dependent boundary condition.
     *
     * \dot{f}_1(x) * (t-t_1(x)) for t >= t_1(x).
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

        const PylithScalar tRel = t - a[i_start];
        if (tRel > 0.0) {
            value[0] = a[i_rate]*tRel;
        } else {
            value[0] = 0.0;
        } // if/else
    } // rate_scalar

    /** Vector rate term for time-dependent boundary condition.
     *
     * \dot{f}_1(x) * (t-t_1(x)) for t >= t_1(x).
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

        const PylithScalar tRel = t - a[i_start];
        if (tRel > 0.0) {
            for (PylithInt i = 0; i < dim; ++i) {
                value[i] = a[i_rate+i]*tRel;
            } // for
        } else {
            for (PylithInt i = 0; i < dim; ++i) {
                value[i] = 0.0;
            } // for
        } // if/else
    } // rate_vector

    /** Scalar time history term for time-dependent boundary condition.
     *
     * f_2(x) * a(t-t_2(x)) for t >= t_2(x).
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
                            const PylithScalar n[],
                            const PylithInt numConstants,
                            const PylithScalar constants[],
                            PylithScalar value[]) {
        const PylithInt _numA = 3;
        assert(_numA == numA);
        assert(aOff);
        assert(a);
        assert(value);

        const PylithInt i_amplitude = aOff[0];
        const PylithInt i_value = aOff[2];

        value[0] = a[i_amplitude]*a[i_value];
    } // timeHistory_scalar

    /** Vector time history term for time-dependent boundary condition.
     *
     * f_2(x) * a(t-t_2(x)) for t >= t_2(x).
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
                            const PylithScalar n[],
                            const PylithInt numConstants,
                            const PylithScalar constants[],
                            PylithScalar value[]) {
        const PylithInt _numA = 3;
        assert(_numA == numA);
        assert(aOff);
        assert(a);
        assert(value);

        const PylithInt i_amplitude = aOff[0];
        const PylithInt i_value = aOff[2];

        for (PylithInt i = 0; i < dim; ++i) {
            value[i] = a[i_amplitude+i]*a[i_value];
        } // for
    } // timeHistory_vector

    /** Compute boundary condition scalar value using initial and rate terms.
     *
     * f_0(x) +
     * \dot{f}_1(x) * H(t-t_1(x)) +
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
                            const PylithScalar n[],
                            const PylithInt numConstants,
                            const PylithScalar constants[],
                            PylithScalar value[]) {
        const PylithInt _numA = 3;
        assert(_numA == numA);
        assert(aOff);
        assert(a);
        assert(value);

        const PylithInt i_initial = 0;
        const PylithInt i_rate = 1;
        const PylithInt i_start = 2;

        const PylithInt numAInitial = 1;
        const PylithInt aOffInitial[1] = { aOff[i_initial] };
        const PylithInt* aOffInitial_x = NULL;

        const PylithInt numARate = 2;
        const PylithInt aOffRate[2] = { aOff[i_rate], aOff[i_start] };
        const PylithInt* aOffRate_x = NULL;

        PylithScalar valueTmp = 0.0;

        initial_scalar(dim, numS, numAInitial, sOff, sOff_x, s, s_t, s_x, aOffInitial, aOffInitial_x, a, a_t, a_x,
                       t, x, n, numConstants, constants, &valueTmp);
        value[0] = valueTmp;

        rate_scalar(dim, numS, numARate, sOff, sOff_x, s, s_t, s_x, aOffRate, aOffRate_x, a, a_t, a_x,
                    t, x, n, numConstants, constants, &valueTmp);
        value[0] += valueTmp;
    } // initialRate_scalar

    /** Compute boundary condition vector value using initial and rate terms.
     *
     * f_0(x) +
     * \dot{f}_1(x) * H(t-t_1(x)) +
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
                            const PylithScalar n[],
                            const PylithInt numConstants,
                            const PylithScalar constants[],
                            PylithScalar value[]) {
        const PylithInt _numA = 3;
        assert(_numA == numA);
        assert(aOff);
        assert(a);
        assert(value);

        const PylithInt i_initial = 0;
        const PylithInt i_rate = 1;
        const PylithInt i_start = 2;

        const PylithInt numAInitial = 1;
        const PylithInt aOffInitial[1] = { aOff[i_initial] };
        const PylithInt* aOffInitial_x = NULL;

        const PylithInt numARate = 2;
        const PylithInt aOffRate[2] = { aOff[i_rate], aOff[i_start] };
        const PylithInt* aOffRate_x = NULL;

        PylithScalar valueTmp[3] = { 0.0, 0.0, 0.0 };

        initial_vector(dim, numS, numAInitial, sOff, sOff_x, s, s_t, s_x, aOffInitial, aOffInitial_x, a, a_t, a_x,
                       t, x, n, numConstants, constants, valueTmp);
        for (PylithInt i = 0; i < dim; ++i) {
            value[i] = valueTmp[i];
        } // for

        rate_vector(dim, numS, numARate, sOff, sOff_x, s, s_t, s_x, aOffRate, aOffRate_x, a, a_t, a_x,
                    t, x, n, numConstants, constants, valueTmp);
        for (PylithInt i = 0; i < dim; ++i) {
            value[i] += valueTmp[i];
        } // for
    } // initialRate_vector

    /** Compute boundary condition scalar value using initial and time history terms.
     *
     * f_0(x) +
     * \dot{f}_1(x) * H(t-t_1(x)) +
     * f_2(x) * a(t-t_2(x)) * H(t-t_2(s).
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
                                   const PylithScalar n[],
                                   const PylithInt numConstants,
                                   const PylithScalar constants[],
                                   PylithScalar value[]) {
        const PylithInt _numA = 4;
        assert(_numA == numA);
        assert(aOff);
        assert(a);
        assert(value);

        const PylithInt i_initial = 0;
        const PylithInt i_thAmp = 1;
        const PylithInt i_thStart = 2;
        const PylithInt i_thValue = 3;

        const PylithInt numAInitial = 1;
        const PylithInt aOffInitial[1] = { aOff[i_initial] };
        const PylithInt* aOffInitial_x = NULL;

        const PylithInt numATimeHistory = 3;
        const PylithInt aOffTimeHistory[3] = { aOff[i_thAmp], aOff[i_thStart], aOff[i_thValue] };
        const PylithInt* aOffTimeHistory_x = NULL;

        PylithScalar valueTmp = 0.0;

        initial_scalar(dim, numS, numAInitial, sOff, sOff_x, s, s_t, s_x, aOffInitial, aOffInitial_x, a, a_t, a_x,
                       t, x, n, numConstants, constants, &valueTmp);
        value[0] = valueTmp;

        timeHistory_scalar(dim, numS, numATimeHistory, sOff, sOff_x, s, s_t, s_x, aOffTimeHistory, aOffTimeHistory_x, a, a_t, a_x,
                           t, x, n, numConstants, constants, &valueTmp);
        value[0] += valueTmp;
    } // initialTimeHistory_scalar

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
                                   const PylithScalar n[],
                                   const PylithInt numConstants,
                                   const PylithScalar constants[],
                                   PylithScalar value[]) {
        const PylithInt _numA = 4;
        assert(_numA == numA);
        assert(aOff);
        assert(a);
        assert(value);

        const PylithInt i_initial = 0;
        const PylithInt i_thAmp = 1;
        const PylithInt i_thStart = 2;
        const PylithInt i_thValue = 3;

        const PylithInt numAInitial = 1;
        const PylithInt aOffInitial[1] = { aOff[i_initial] };
        const PylithInt* aOffInitial_x = NULL;

        const PylithInt numATimeHistory = 3;
        const PylithInt aOffTimeHistory[3] = { aOff[i_thAmp], aOff[i_thStart], aOff[i_thValue] };
        const PylithInt* aOffTimeHistory_x = NULL;

        PylithScalar valueTmp[3] = { 0.0, 0.0, 0.0 };

        initial_vector(dim, numS, numAInitial, sOff, sOff_x, s, s_t, s_x, aOffInitial, aOffInitial_x, a, a_t, a_x,
                       t, x, n, numConstants, constants, valueTmp);
        for (PylithInt i = 0; i < dim; ++i) {
            value[i] = valueTmp[i];
        } // for

        timeHistory_vector(dim, numS, numATimeHistory, sOff, sOff_x, s, s_t, s_x, aOffTimeHistory, aOffTimeHistory_x, a, a_t, a_x,
                           t, x, n, numConstants, constants, valueTmp);
        for (PylithInt i = 0; i < dim; ++i) {
            value[i] += valueTmp[i];
        } // for
    } // initialTimeHistory_vector

    /** Compute boundary condition scalar value using rate and time history terms.
     *
     * \dot{f}_1(x) * H(t-t_1(x)) +
     * f_2(x) * a(t-t_2(x)) * H(t-t_2(s).
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
                                const PylithScalar n[],
                                const PylithInt numConstants,
                                const PylithScalar constants[],
                                PylithScalar value[]) {
        const PylithInt _numA = 5;
        assert(_numA == numA);
        assert(aOff);
        assert(a);
        assert(value);

        const PylithInt i_rateAmp = 0;
        const PylithInt i_rateStart = 0;
        const PylithInt i_thAmp = 1;
        const PylithInt i_thStart = 2;
        const PylithInt i_thValue = 3;

        const PylithInt numARate = 1;
        const PylithInt aOffRate[2] = { aOff[i_rateAmp], aOff[i_rateStart] };
        const PylithInt* aOffRate_x = NULL;

        const PylithInt numATimeHistory = 3;
        const PylithInt aOffTimeHistory[3] = { aOff[i_thAmp], aOff[i_thStart], aOff[i_thValue] };
        const PylithInt* aOffTimeHistory_x = NULL;

        PylithScalar valueTmp = 0.0;

        rate_scalar(dim, numS, numARate, sOff, sOff_x, s, s_t, s_x, aOffRate, aOffRate_x, a, a_t, a_x,
                    t, x, n, numConstants, constants, &valueTmp);
        value[0] = valueTmp;

        timeHistory_scalar(dim, numS, numATimeHistory, sOff, sOff_x, s, s_t, s_x, aOffTimeHistory, aOffTimeHistory_x, a, a_t, a_x,
                           t, x, n, numConstants, constants, &valueTmp);
        value[0] += valueTmp;
    } // rateTimeHistory_scalar

    /** Compute boundary condition vector value using rate and time history terms.
     *
     * \dot{f}_1(x) * H(t-t_1(x)) +
     * f_2(x) * a(t-t_2(x)) * H(t-t_2(s).
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
                                const PylithScalar n[],
                                const PylithInt numConstants,
                                const PylithScalar constants[],
                                PylithScalar value[]) {
        const PylithInt _numA = 5;
        assert(_numA == numA);
        assert(aOff);
        assert(a);
        assert(value);

        const PylithInt i_rateAmp = 0;
        const PylithInt i_rateStart = 0;
        const PylithInt i_thAmp = 1;
        const PylithInt i_thStart = 2;
        const PylithInt i_thValue = 3;

        const PylithInt numARate = 1;
        const PylithInt aOffRate[2] = { aOff[i_rateAmp], aOff[i_rateStart] };
        const PylithInt* aOffRate_x = NULL;

        const PylithInt numATimeHistory = 3;
        const PylithInt aOffTimeHistory[3] = { aOff[i_thAmp], aOff[i_thStart], aOff[i_thValue] };
        const PylithInt* aOffTimeHistory_x = NULL;

        PylithScalar valueTmp[3] = { 0.0, 0.0, 0.0 };

        rate_vector(dim, numS, numARate, sOff, sOff_x, s, s_t, s_x, aOffRate, aOffRate_x, a, a_t, a_x,
                    t, x, n, numConstants, constants, valueTmp);
        for (PylithInt i = 0; i < dim; ++i) {
            value[i] = valueTmp[i];
        } // for

        timeHistory_vector(dim, numS, numATimeHistory, sOff, sOff_x, s, s_t, s_x,
                           aOffTimeHistory, aOffTimeHistory_x, a, a_t, a_x,
                           t, x, n, numConstants, constants, valueTmp);
        for (PylithInt i = 0; i < dim; ++i) {
            value[i] += valueTmp[i];
        } // for
    } // rateTimeHistory_vector

    /** Compute boundary condition scalar value using initial, rate ,and time history terms.
     *
     * f_0(x) +
     * \dot{f}_1(x) * H(t-t_1(x)) +
     * f_2(x) * a(t-t_2(x)) * H(t-t_2(s).
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
                                       const PylithScalar n[],
                                       const PylithInt numConstants,
                                       const PylithScalar constants[],
                                       PylithScalar value[]) {
        const PylithInt _numA = 6;
        assert(_numA == numA);
        assert(aOff);
        assert(a);
        assert(value);

        const PylithInt i_valueAmp = 0;
        const PylithInt i_rateAmp = 1;
        const PylithInt i_rateStart = 2;
        const PylithInt i_thAmp = 3;
        const PylithInt i_thStart = 4;
        const PylithInt i_thValue = 5;

        const PylithInt numAInitial = 1;
        const PylithInt aOffInitial[1] = { aOff[i_valueAmp] };
        const PylithInt* aOffInitial_x = NULL;

        const PylithInt numARate = 1;
        const PylithInt aOffRate[2] = { aOff[i_rateAmp], aOff[i_rateStart] };
        const PylithInt* aOffRate_x = NULL;

        const PylithInt numATimeHistory = 3;
        const PylithInt aOffTimeHistory[3] = { aOff[i_thAmp], aOff[i_thStart], aOff[i_thValue] };
        const PylithInt* aOffTimeHistory_x = NULL;

        PylithScalar valueTmp = 0.0;

        initial_scalar(dim, numS, numAInitial, sOff, sOff_x, s, s_t, s_x, aOffInitial, aOffInitial_x, a, a_t, a_x,
                       t, x, n, numConstants, constants, &valueTmp);
        value[0] = valueTmp;

        rate_scalar(dim, numS, numARate, sOff, sOff_x, s, s_t, s_x, aOffRate, aOffRate_x, a, a_t, a_x,
                    t, x, n, numConstants, constants, &valueTmp);
        value[0] += valueTmp;

        timeHistory_scalar(dim, numS, numATimeHistory, sOff, sOff_x, s, s_t, s_x,
                           aOffTimeHistory, aOffTimeHistory_x, a, a_t, a_x, t, x, n, numConstants, constants, &valueTmp);
        value[0] += valueTmp;
    } // initialRateTimeHistory_scalar

    /** Compute boundary condition vector value using initial, rate, and time history terms.
     *
     * f_0(x) +
     * \dot{f}_1(x) * H(t-t_1(x)) +
     * f_2(x) * a(t-t_2(x)) * H(t-t_2(s).
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
                                       const PylithScalar n[],
                                       const PylithInt numConstants,
                                       const PylithScalar constants[],
                                       PylithScalar value[]) {
        const PylithInt _numA = 6;
        assert(_numA == numA);
        assert(aOff);
        assert(a);
        assert(value);

        const PylithInt i_valueAmp = 0;
        const PylithInt i_rateAmp = 1;
        const PylithInt i_rateStart = 2;
        const PylithInt i_thAmp = 3;
        const PylithInt i_thStart = 4;
        const PylithInt i_thValue = 5;

        const PylithInt numAInitial = 1;
        const PylithInt aOffInitial[1] = { aOff[i_valueAmp] };
        const PylithInt* aOffInitial_x = NULL;

        const PylithInt numARate = 1;
        const PylithInt aOffRate[2] = { aOff[i_rateAmp], aOff[i_rateStart] };
        const PylithInt* aOffRate_x = NULL;

        const PylithInt numATimeHistory = 3;
        const PylithInt aOffTimeHistory[3] = { aOff[i_thAmp], aOff[i_thStart], aOff[i_thValue] };
        const PylithInt* aOffTimeHistory_x = NULL;

        PylithScalar valueTmp[3] = { 0.0, 0.0, 0.0 };

        initial_vector(dim, numS, numAInitial, sOff, sOff_x, s, s_t, s_x, aOffInitial, aOffInitial_x, a, a_t, a_x,
                       t, x, n, numConstants, constants, valueTmp);
        for (PylithInt i = 0; i < dim; ++i) {
            value[i] = valueTmp[i];
        } // for

        rate_vector(dim, numS, numARate, sOff, sOff_x, s, s_t, s_x, aOffRate, aOffRate_x, a, a_t, a_x,
                    t, x, n, numConstants, constants, valueTmp);
        for (PylithInt i = 0; i < dim; ++i) {
            value[i] += valueTmp[i];
        } // for

        timeHistory_vector(dim, numS, numATimeHistory, sOff, sOff_x, s, s_t, s_x,
                           aOffTimeHistory, aOffTimeHistory_x, a, a_t, a_x, t, x, n, numConstants, constants, valueTmp);
        for (PylithInt i = 0; i < dim; ++i) {
            value[i] += valueTmp[i];
        } // for
    } // initialRateTimeHistory_vector

}; // TimeDependentFn

#endif // pylith_fekernels_timedependentfn_hh

// End of file
