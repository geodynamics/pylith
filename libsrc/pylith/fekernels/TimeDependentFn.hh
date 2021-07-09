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
 * Copyright (c) 2010-2021 University of California, Davis
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
    static
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
                        PylithScalar value[]);

    /** Vector initial value term for time-dependent boundary condition.
     *
     * f_0(x)
     */
    static
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
                        PylithScalar value[]);

    /** Scalar rate term for time-dependent boundary condition.
     *
     * \dot{f}_1(x) * (t-t_1(x)) for t >= t_1(x).
     */
    static
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
                     PylithScalar value[]);

    /** Vector rate term for time-dependent boundary condition.
     *
     * \dot{f}_1(x) * (t-t_1(x)) for t >= t_1(x).
     */
    static
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
                     PylithScalar value[]);

    /** Scalar time history term for time-dependent boundary condition.
     *
     * f_2(x) * a(t-t_2(x)) for t >= t_2(x).
     */
    static
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
                            PylithScalar value[]);

    /** Vector time history term for time-dependent boundary condition.
     *
     * f_2(x) * a(t-t_2(x)) for t >= t_2(x).
     */
    static
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
                            PylithScalar value[]);

    /** Compute boundary condition scalar value using initial and rate terms.
     *
     * f_0(x) +
     * \dot{f}_1(x) * H(t-t_1(x)) +
     */
    static
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
                            PylithScalar value[]);

    /** Compute boundary condition vector value using initial and rate terms.
     *
     * f_0(x) +
     * \dot{f}_1(x) * H(t-t_1(x)) +
     */
    static
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
                            PylithScalar value[]);

    /** Compute boundary condition scalar value using initial and time history terms.
     *
     * f_0(x) +
     * \dot{f}_1(x) * H(t-t_1(x)) +
     * f_2(x) * a(t-t_2(x)) * H(t-t_2(s).
     */
    static
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
                                   PylithScalar value[]);

    /** Compute boundary condition vector value using initial and time history terms.
     *
     * f_0(x) +
     * \dot{f}_1(x) * H(t-t_1(x)) +
     * f_2(x) * a(t-t_2(x)) * H(t-t_2(s).
     */
    static
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
                                   PylithScalar value[]);

    /** Compute boundary condition scalar value using rate and time history terms.
     *
     * \dot{f}_1(x) * H(t-t_1(x)) +
     * f_2(x) * a(t-t_2(x)) * H(t-t_2(s).
     */
    static
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
                                PylithScalar value[]);

    /** Compute boundary condition vector value using rate and time history terms.
     *
     * \dot{f}_1(x) * H(t-t_1(x)) +
     * f_2(x) * a(t-t_2(x)) * H(t-t_2(s).
     */
    static
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
                                PylithScalar value[]);

    /** Compute boundary condition scalar value using initial, rate ,and time history terms.
     *
     * f_0(x) +
     * \dot{f}_1(x) * H(t-t_1(x)) +
     * f_2(x) * a(t-t_2(x)) * H(t-t_2(s).
     */
    static
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
                                       PylithScalar value[]);

    /** Compute boundary condition vector value using initial, rate, and time history terms.
     *
     * f_0(x) +
     * \dot{f}_1(x) * H(t-t_1(x)) +
     * f_2(x) * a(t-t_2(x)) * H(t-t_2(s).
     */
    static
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
                                       PylithScalar value[]);

}; // TimeDependentFn

#endif // pylith_fekernels_timedependentfn_hh


// End of file
