/*
 * ================================================================================================
 * This code is part of PyLith, developed through the Computational Infrastructure
 * for Geodynamics (https://github.com/geodynamics/pylith).
 *
 * Copyright (c) 2010-2023, University of California, Davis and the PyLith Development Team.
 * All rights reserved.
 *
 * See https://mit-license.org/ and LICENSE.md and for license information.
 * =================================================================================================
 */

/** @file libsrc/fekernels/AbsorbingDampers.hh
 *
 * Kernels for absorbing dampers boundary conditions.
 *
 * \int_{\Gamma_\tau} \trialvec[u] \vec{\tau}(\vec{x},t) d\Gamma
 *
 * where \tau_{shear} = -rho(\vec{x}) v_s v_{shear}
 * and \tau_{normal} = -rho(\vec{x}) v_p v_{normal}
 *
 * Solution fields: [disp(dim), vel(dim), ...]
 *
 * Auxiliary fields:
 * - 0: density(1)
 * - 1: vp(1)
 * - 2: vs(1)
 */

#if !defined(pylith_fekernels_absorbingdampers_hh)
#define pylith_fekernels_absorbingdampers_hh

// Include directives ---------------------------------------------------
#include "fekernelsfwd.hh" // forward declarations

#include "pylith/utils/types.hh"

#include <cassert> // USES assert()

class pylith::fekernels::AbsorbingDampers {
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

    /** g0 function for absorbing dampers boundary condition.
     *
     * g_0(x)
     */
    static inline
    void g0(const PylithInt dim,
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
            PylithScalar g0[]) {
        assert(2 == dim || 3 == dim);

        const PylithInt _numA = 3;
        assert(_numA == numA);
        assert(aOff);
        assert(a);
        const PylithInt i_density = aOff[0];
        const PylithInt i_vp = aOff[1];
        const PylithInt i_vs = aOff[2];

        const PylithInt _numS = 2;
        assert(sOff);
        assert(s);
        assert(numS >= _numS);
        const PylithInt i_vel = sOff[1];

        const PylithScalar density = a[i_density];
        const PylithScalar vp = a[i_vp];
        const PylithScalar vs = a[i_vs];

        PylithScalar velN[3];
        PylithScalar velT[3];
        PylithScalar velNMag = 0;
        for (PylithInt i = 0; i < dim; ++i) {
            velNMag += s[i_vel+i] * n[i];
        } // for
        for (PylithInt i = 0; i < dim; ++i) {
            velN[i] = velNMag * n[i];
        } // for
        for (PylithInt i = 0; i < dim; ++i) {
            velT[i] = s[i_vel+i] - velN[i];
        } // for

        for (PylithInt i = 0; i < dim; ++i) {
            g0[i] -= density * (vs * velT[i] + vp * velN[i]);
        } // for
    } // g0

}; // AbsorbingDampers

#endif // pylith_fekernels_absorbingdampers_hh

// End of file
