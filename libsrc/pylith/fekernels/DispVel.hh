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

/*
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

// Include directives ---------------------------------------------------
#include "pylith/fekernels/fekernelsfwd.hh" // forward declarations

#include "pylith/utils/types.hh"

#include <cassert> // USES assert()

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
    static inline
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
             PylithScalar f0[]) {
        assert(sOff);
        assert(s);
        assert(s_t);
        assert(f0);

        const PylithInt _numS = 2;
        assert(_numS == numS);

        const PylithInt i_disp = 0;
        const PylithScalar* disp_t = &s_t[sOff[i_disp]];

        for (PylithInt i = 0; i < dim; ++i) {
            f0[i] += disp_t[i];
        } // for
    } // f0u

    /** f0 function for velocity equation: f0u = \dot{v}.
     *
     * Solution fields: [disp(dim), vel(dim)]
     */
    static inline
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
             PylithScalar f0[]) {
        assert(sOff);
        assert(s);
        assert(s_t);
        assert(f0);

        const PylithInt _numS = 2;
        assert(_numS == numS);

        const PylithInt i_vel = 1;
        const PylithScalar* vel_t = &s_t[sOff[i_vel]];

        for (PylithInt i = 0; i < dim; ++i) {
            f0[i] += vel_t[i];
        } // for
    } // f0v

    /** g0 function for displacement equation: g0u = v.
     *
     * Solution fields: [disp(dim), vel(dim)]
     */
    static inline
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
             PylithScalar g0[]) {
        assert(sOff);
        assert(s);
        assert(g0);

        const PylithInt _numS = 2;
        assert(_numS == numS);

        const PylithInt i_vel = 1;
        const PylithScalar* vel = &s[sOff[i_vel]];

        for (PylithInt i = 0; i < dim; ++i) {
            g0[i] += vel[i];
        } // for
    } // g0u

    /** Jf0 function for displacement equation with zero values on diagonal.
     *
     * This is associated with the elasticity equation without intertia.
     *
     * Solution fields: [...]
     */
    static inline
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
                    PylithScalar Jf0[]) {
        // No work to do for zero values.
    } // Jf0uu_zero

    /** Jf0 function for displacement equation: Jf0uu = s_tshift.
     *
     * Solution fields: [disp(dim), vel(dim)]
     */
    static inline
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
                       PylithScalar Jf0[]) {
        const PylithInt _numS = 2;
        assert(_numS == numS);
        assert(s_tshift > 0);

        for (PylithInt i = 0; i < dim; ++i) {
            Jf0[i*dim+i] += s_tshift;
        } // for
    } // Jf0uu_stshift

    /** Jg0 function for displacement equation: 1.0.
     *
     * Solution fields: [disp(dim), vel(dim)]
     */
    static inline
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
               PylithScalar Jg0[]) {
        const PylithInt _numS = 2;
        assert(_numS == numS);

        for (PylithInt i = 0; i < dim; ++i) {
            Jg0[i*dim+i] += 1.0;
        } // for
    } // Jg0uv

}; // DispVel

/* End of file */
