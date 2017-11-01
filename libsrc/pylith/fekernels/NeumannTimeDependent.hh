/* -*- C -*-
 *
 * ----------------------------------------------------------------------
 *
 * Brad T. Aagaard, U.S. Geological Survey
 * Charles A. Williams, GNS Science
 * Matthew G. Knepley, University of Chicago
 *
 * This code was developed as part of the Computational Infrastructure
 * for Geodynamics (http:*geodynamics.org).
 *
 * Copyright (c) 2010-2017 University of California, Davis
 *
 * See COPYING for license information.
 *
 * ----------------------------------------------------------------------
 */

/** @file libsrc/fekernels/NeumannTimeDependent.h
 *
 * Kernels for computing value from parameters for Neumann time dependent boundary conditions.
 *
 * \int_{\Gamma_\tau} \trialvec[u] \vec{\tau}(\vec{x},t) d\Gamma
 */

#if !defined(pylith_fekernels_neumanntimedependent_hh)
#define pylith_fekernels_neumanntimedependent_hh

// Include directives ---------------------------------------------------
#include "fekernelsfwd.hh" // forward declarations

#include "pylith/utils/types.hh"

class pylith::fekernels::NeumannTimeDependent {

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
     * @param[in] n Unit vector normal to boundary.
     * @param[in] numConstants Number of registered constants.
     * @param[in] constants Array of registered constants.
     * @param[out] g0 [dim].
     */

    /** Scalar initial value term for time-dependent boundary condition.
     *
     * f_0(x)
     */
    static
    void g0_initial_scalar(const PylithInt dim,
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
                           PylithScalar g0[]);

    /** Vector initial value term for time-dependent boundary condition.
     *
     * f_0(x)
     */
    static
    void g0_initial_vector(const PylithInt dim,
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
                           PylithScalar g0[]);

    /** Scalar rate term for time-dependent boundary condition.
     *
     * \dot{f}_1(x) * (t-t_1(x)) for t >= t_1(x).
     */
    static
    void g0_rate_scalar(const PylithInt dim,
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
                        PylithScalar g0[]);

    /** Vector rate term for time-dependent boundary condition.
     *
     * \dot{f}_1(x) * (t-t_1(x)) for t >= t_1(x).
     */
    static
    void g0_rate_vector(const PylithInt dim,
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
                        PylithScalar g0[]);

    /** Scalar time history term for time-dependent boundary condition.
     *
     * f_2(x) * a(t-t_2(x)) for t >= t_2(x).
     */
    static
    void g0_timeHistory_scalar(const PylithInt dim,
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
                               PylithScalar g0[]);

    /** Vector time history term for time-dependent boundary condition.
     *
     * f_2(x) * a(t-t_2(x)) for t >= t_2(x).
     */
    static
    void g0_timeHistory_vector(const PylithInt dim,
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
                               PylithScalar g0[]);

    /** Compute boundary condition scalar value using initial and rate terms.
     *
     * f_0(x) +
     * \dot{f}_1(x) * H(t-t_1(x)) +
     */
    static
    void g0_initialRate_scalar(const PylithInt dim,
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
                               PylithScalar g0[]);

    /** Compute boundary condition vector value using initial and rate terms.
     *
     * f_0(x) +
     * \dot{f}_1(x) * H(t-t_1(x)) +
     */
    static
    void g0_initialRate_vector(const PylithInt dim,
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
                               PylithScalar g0[]);

    /** Compute boundary condition scalar value using initial and time history terms.
     *
     * f_0(x) +
     * \dot{f}_1(x) * H(t-t_1(x)) +
     * f_2(x) * a(t-t_2(x)) * H(t-t_2(s).
     */
    static
    void g0_initialTimeHistory_scalar(const PylithInt dim,
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
                                      PylithScalar g0[]);

    /** Compute boundary condition vector value using initial and time history terms.
     *
     * f_0(x) +
     * \dot{f}_1(x) * H(t-t_1(x)) +
     * f_2(x) * a(t-t_2(x)) * H(t-t_2(s).
     */
    static
    void g0_initialTimeHistory_vector(const PylithInt dim,
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
                                      PylithScalar g0[]);

    /** Compute boundary condition scalar value using rate and time history terms.
     *
     * \dot{f}_1(x) * H(t-t_1(x)) +
     * f_2(x) * a(t-t_2(x)) * H(t-t_2(s).
     */
    static
    void g0_rateTimeHistory_scalar(const PylithInt dim,
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
                                   PylithScalar g0[]);

    /** Compute boundary condition vector value using rate and time history terms.
     *
     * \dot{f}_1(x) * H(t-t_1(x)) +
     * f_2(x) * a(t-t_2(x)) * H(t-t_2(s).
     */
    static
    void g0_rateTimeHistory_vector(const PylithInt dim,
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
                                   PylithScalar g0[]);

    /** Compute boundary condition scalar value using initial, rate ,and time history terms.
     *
     * f_0(x) +
     * \dot{f}_1(x) * H(t-t_1(x)) +
     * f_2(x) * a(t-t_2(x)) * H(t-t_2(s).
     */
    static
    void g0_initialRateTimeHistory_scalar(const PylithInt dim,
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
                                          PylithScalar g0[]);

    /** Compute boundary condition vector value using initial, rate, and time history terms.
     *
     * f_0(x) +
     * \dot{f}_1(x) * H(t-t_1(x)) +
     * f_2(x) * a(t-t_2(x)) * H(t-t_2(s).
     */
    static
    void g0_initialRateTimeHistory_vector(const PylithInt dim,
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
                                          PylithScalar g0[]);

}; // NeumannTimeDependent

#endif // pylith_fekernels_neumanntimedependent_hh


// End of file
