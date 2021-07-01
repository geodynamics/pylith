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

/** @file libsrc/fekernels/DispVel.hh
 *
 * Kernels for time evolution equation with displacement and velocity
 * solution fields.
 *
 * Solution fields: [disp(dim), vel(dim, optional)]
 *
 * Auxiliary fields: [...] (not used)
 *
 * \int_V \vec{\phi}_v \cdot \left( \frac{\partial \vec{u}(t)}{\partial t} \right) \, dV =
 *   \int_V \vec{\phi}_v \cdot \vec{v}(t) \, dV.
 *
 * LHS Residual
 *
 * f0_DispVel: \vec{f0} = \frac{\partial \vec{u}(t)}{\partial t}
 *
 * RHS Residual
 *
 * g0_DispVel: \vec{g0} = \vec{v}(t)
 *
 * LHS Jacobian
 *
 * Jf0_veldisp_DispVelImplicit: s_tshift
 *
 * Jf0_veldisp_DispVelExplicit: 0
 *
 * RHS Jacobian
 *
 * Jg0_velvel_DispVel: +1.0
 *
 * ======================================================================
 */

#if !defined(pylith_fekernels_DispVel_hh)
#define pylith_fekernels_DispVel_hh

// Include directives ---------------------------------------------------
#include "fekernelsfwd.hh" // forward declarations

#include "pylith/utils/types.hh"

class pylith::fekernels::DispVel {

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


    /** f0 function for displacement equation: f0u = \dot{u}.
     *
     * Solution fields: [disp(dim), vel(dim)]
     */
    static
    void f0u(const PylithInt dim,
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
             PylithScalar f0[]);


    /** f0 function for velocity equation: f0u = \dot{v}.
     *
     * Solution fields: [disp(dim), vel(dim)]
     */
    static
    void f0v(const PylithInt dim,
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
             PylithScalar f0[]);


    /** g0 function for displacement equation: g0u = v.
     *
     * Solution fields: [disp(dim), vel(dim)]
     */
    static
    void g0u(const PylithInt dim,
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
             PylithScalar g0[]);


    /** Jf0 function for displacement equation with zero values on diagonal.
     *
     * This is associated with the elasticity equation without intertia.
     *
     * Solution fields: [...]
     */
    static
    void Jf0uu_zero(const PylithInt dim,
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
                    const PylithInt numConstants,
                    const PylithScalar constants[],
                    PylithScalar Jf0[]);

    /** Jf0 function for displacement equation: Jf0uu = s_tshift.
     *
     * Solution fields: [disp(dim), vel(dim)]
     */
    static
    void Jf0uu_stshift(const PylithInt dim,
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
                       const PylithInt numConstants,
                       const PylithScalar constants[],
                       PylithScalar Jf0[]);

    /** Jg0 function for displacement equation: 1.0.
     *
     * Solution fields: [disp(dim), vel(dim)]
     */
    static
    void Jg0uv(const PylithInt dim,
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
               const PylithInt numConstants,
               const PylithScalar constants[],
               PylithScalar Jg0[]);


}; // DispVel

#endif /* pylith_fekernels_DispVel_hh */


/* End of file */
