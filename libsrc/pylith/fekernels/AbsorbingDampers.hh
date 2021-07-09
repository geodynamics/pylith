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
    static
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
            PylithScalar g0[]);

}; // AbsorbingDampers

#endif // pylith_fekernels_absorbingdampers_hh


// End of file
