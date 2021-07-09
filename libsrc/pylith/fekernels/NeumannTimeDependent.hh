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

private:
    // PRIVATE MEMBERS ///////////////////////////////////////////////////////

    /** Compute tangential directions from reference direction (first and second choice) and normal direction in 3-D.
     *
     * We use the first reference direction if it is not colinear with the normal direction; otherwise, we use the second
     * reference direction.
     *
     * + tangential direction 1 is refDir x normal.
     * + tangential direction 2 is normal x tangential direction 1.
     *
     * @param[in] dim Dimension (must be 3)
     * @param[in] refDir1 Reference direction (first choice) [3].
     * @param[in] refDir2 Reference direction (second choice) [3].
     * @param[in] normDir Normal direction [3].
     * @param[out] tanDir1 First tangential direction [3].
     * @param[out] tanDir2 Second tangential direction [3].
     */
    static
    void _tangential_directions(const PylithInt dim,
                                const PylithScalar refDir1[],
                                const PylithScalar refDir2[],
                                const PylithScalar normDir[],
                                PylithScalar tanDir1[],
                                PylithScalar tanDir2[]);

    // PUBLIC MEMBERS ///////////////////////////////////////////////////////
}; // NeumannTimeDependent

#endif // pylith_fekernels_neumanntimedependent_hh


// End of file
