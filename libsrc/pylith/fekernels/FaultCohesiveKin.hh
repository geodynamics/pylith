/* -*- C++ -*-
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

/** @file libsrc/fekernels/FaultCohesiveKin.hh
 *
 * Kernels for faults with prescribed slip.
 *
 * Solution fields: [disp(dim), vel(dim, optional), lagrange(dim)]
 *
 * Auxiliary fields: [...] (not used)
 *
 * LHS Residual: no contribution
 *
 * RHS Residual
 *
 *  - g0u^+ = -\lambda
 *  - g0u^- = +\lambda
 *  - g0\lambda = d - u^+ + u^-
 *
 * LHS Jacobian: no contribution
 *
 * RHS Jacobian
 *
 *  - Jg0^{u \lambda} = -1 for u^+, +1 for u^-
 *  - Jg0^{\lambda u} = -1 for u^+, +1 for u^-
 *
 * ======================================================================
 */

#if !defined(pylith_fekernels_faultcohesivekin_hh)
#define pylith_fekernels_faultcohesivekin_hh

// Include directives ---------------------------------------------------
#include "fekernelsfwd.hh" // forward declarations

#include "pylith/utils/types.hh"

class pylith::fekernels::FaultCohesiveKin {
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
     * @param[in] numConstants Number of registered constants.
     * @param[in] constants Array of registered constants.
     * @param[out] f0 [dim].
     */

    /** f0 function for elasticity equation: f0u = +\lambda (neg side).
     *
     * Solution fields: [disp(dim), ..., lagrange(dim)]
     */
    static
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
                 PylithScalar f0[]);

    /** f0 function for elasticity equation: f0u = -\lambda (pos side).
     *
     * Solution fields: [disp(dim), ..., lagrange(dim)]
     */
    static
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
                 PylithScalar f0[]);

    /** f0 function for slip constraint equation: f0\lambda = (u^+ - u^-) - d
     *
     * Solution fields: [disp(dim), ..., lagrange(dim)]
     */
    static
    void f0l_u(const PylithInt dim,
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
               PylithScalar f0[]);

    /** f0 function for slip acceleration constraint equation: f0\lambda = (\dot{v}^+ - \dot{v}^-) - \ddot{d}
     *
     * Solution fields: [disp(dim), vel(dim), ..., lagrange(dim)]
     */
    static
    void f0l_a(const PylithInt dim,
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
               PylithScalar f0[]);

    /** Jf0 function for displacement equation: -\lambda (neg side).
     */
    static
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
                   PylithScalar Jf0[]);

    /** Jf0 function for displacement equation: +\lambda (pos side).
     */
    static
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
                   PylithScalar Jf0[]);

    /** Jf0 function for slip constraint equation: +\lambda (pos side), -\lambda (neg side).
     *
     * Solution fields: [disp(dim), ..., lagrange(dim)]
     */
    static
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
               PylithScalar Jf0[]);

}; // FaultCohesiveKin

#endif // pylith_fekernels_faultcohesivekin_hh

/* End of file */
