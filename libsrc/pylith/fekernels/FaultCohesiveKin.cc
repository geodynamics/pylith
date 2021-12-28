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

#include <portinfo>

#include "pylith/fekernels/FaultCohesiveKin.hh"

#include "pylith/fekernels/BoundaryDirections.hh" // USES tangential_directions()

#include <cassert> // USES assert()

/* ======================================================================
 * Kernels for prescribed fault slip.
 *
 * Solution fields = [disp(dim), vel(dim) [if dynamic], lagrange_multiplier]
 * Auxiliary fields = [slip (dim)]
 *
 * ======================================================================
 */

// ----------------------------------------------------------------------
// f0 function for elasticity equation: f0u = +\lambda (neg side).
void
pylith::fekernels::FaultCohesiveKin::f0u_neg(const PylithInt dim,
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
                                             PylithScalar f0[]) {
    assert(sOff);
    assert(s);
    assert(f0);

    assert(numS >= 2);

    const PylithInt spaceDim = dim + 1; // :KLUDGE: dim passed in is spaceDim-1

    const PylithInt sOffLagrange = sOff[numS-1];
    const PylithScalar* lagrange = &s[sOffLagrange];

    for (PylithInt i = 0; i < spaceDim; ++i) {
        f0[i] += -lagrange[i];
    } // for
} // f0u_neg


// ----------------------------------------------------------------------
// f0 function for elasticity equation: f0u = -\lambda (pos side).
void
pylith::fekernels::FaultCohesiveKin::f0u_pos(const PylithInt dim,
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
                                             PylithScalar f0[]) {
    assert(sOff);
    assert(s);
    assert(f0);

    assert(numS >= 2);

    const PylithInt spaceDim = dim + 1; // :KLUDGE: dim passed in is spaceDim-1

    const PylithInt sOffLagrange = sOff[numS-1];
    const PylithScalar* lagrange = &s[sOffLagrange];

    for (PylithInt i = 0; i < spaceDim; ++i) {
        // f0[fOffN+i] += -lagrange[i];
        f0[i] += +lagrange[i];
    } // for
} // f0u_pos


#include <iostream>
// ----------------------------------------------------------------------
// f0 function for slip constraint equation: f0\lambda = (u^+ - u^-) - d
void
pylith::fekernels::FaultCohesiveKin::f0l_u(const PylithInt dim,
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
                                           PylithScalar f0[]) {
    assert(sOff);
    assert(aOff);
    assert(s);
    assert(a);
    assert(f0);

    assert(numS >= 2);
    assert(numA >= 1);

    const PylithInt spaceDim = dim + 1; // :KLUDGE: dim passed in is spaceDim-1
    const PylithInt i_slip = 0;
    const PylithInt i_disp = 0;

    const PylithScalar* slip = &a[aOff[i_slip]];

    const PylithInt sOffDispN = sOff[i_disp];
    const PylithInt sOffDispP = sOffDispN+spaceDim;

    const PylithScalar* dispN = &s[sOffDispN];
    const PylithScalar* dispP = &s[sOffDispP];

    switch (spaceDim) {
    case 2: {
        const PylithInt _spaceDim = 2;
        const PylithScalar tanDir[2] = {-n[1], n[0] };
        for (PylithInt i = 0; i < _spaceDim; ++i) {
            const PylithScalar slipXY = n[i]*slip[0] + tanDir[i]*slip[1];
            f0[i] += dispP[i] - dispN[i] - slipXY;
        } // for
        std::cout<<"t="<<t<<", slip="<<slip[1]<<", relDisp="<<dispP[1]-dispN[1]<<std::endl;
        break;
    } // case 2
    case 3: {
        const PylithInt _spaceDim = 3;
        const PylithScalar* refDir1 = &constants[0];
        const PylithScalar* refDir2 = &constants[3];
        PylithScalar tanDir1[3], tanDir2[3];
        BoundaryDirections::tangential_directions(tanDir1, tanDir2, refDir1, refDir2, n);

        for (PylithInt i = 0; i < _spaceDim; ++i) {
            const PylithScalar slipXYZ = n[i]*slip[0] + tanDir1[i]*slip[1] + tanDir2[i]*slip[2];
            f0[i] += dispP[i] - dispN[i] - slipXYZ;
        } // for
        break;
    } // case 3
    default:
        assert(0);
    } // switch
} // f0l_u


// ----------------------------------------------------------------------
// f0 function for slip acceleration constraint equation: f0\lambda = (\dot{v}^+ - \dot{v}^-) - \ddot{d}
void
pylith::fekernels::FaultCohesiveKin::f0l_a(const PylithInt dim,
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
                                           PylithScalar f0[]) {
    assert(sOff);
    assert(aOff);
    assert(s);
    assert(a);
    assert(f0);

    assert(numS >= 3);
    assert(numA >= 1);

    const PylithInt spaceDim = dim + 1; // :KLUDGE: dim passed in is spaceDim-1
    const PylithInt i_slipAcc = numA-1;

    const PylithScalar* slipAcc = &a[aOff[i_slipAcc]];

    switch (spaceDim) {
    case 2: {
        const PylithInt _spaceDim = 2;
        const PylithScalar tanDir[2] = {-n[1], n[0] };
        for (PylithInt i = 0; i < _spaceDim; ++i) {
            const PylithScalar slipAccXY = n[i]*slipAcc[0] + tanDir[i]*slipAcc[1];
            f0[i] += slipAccXY;
        } // for
        break;
    } // case 2
    case 3: {
        const PylithInt _spaceDim = 3;
        const PylithScalar* refDir1 = &constants[0];
        const PylithScalar* refDir2 = &constants[3];
        PylithScalar tanDir1[3], tanDir2[3];
        BoundaryDirections::tangential_directions(tanDir1, tanDir2, refDir1, refDir2, n);

        for (PylithInt i = 0; i < _spaceDim; ++i) {
            const PylithScalar slipAccXYZ = n[i]*slipAcc[0] + tanDir1[i]*slipAcc[1] + tanDir2[i]*slipAcc[2];
            f0[i] += slipAccXYZ;
        } // for
        break;
    } // case 3
    default:
        assert(0);
    } // switch
} // f0l_a


// ----------------------------------------------------------------------
/* Jf0 function for integration of the displacement equation (neg side).
 *
 * Solution fields = [disp(dim), ..., lagrange(dim)]
 * Auxiliary fields = None
 */
void
pylith::fekernels::FaultCohesiveKin::Jf0ul_neg(const PylithInt dim,
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
                                               PylithScalar Jf0[]) {
    assert(numS >= 2);
    assert(Jf0);
    assert(sOff);
    assert(n);

    const PylithInt spaceDim = dim + 1; // :KLUDGE: dim passed in is spaceDim-1

    const PylithInt ncols = spaceDim;

    for (PylithInt i = 0; i < spaceDim; ++i) {
        Jf0[i*ncols+i] += -1.0;
    } // for
} // Jg0ul_neg


// ----------------------------------------------------------------------
/* Jf0 function for integration of the displacement equation (pos side).
 *
 * Solution fields = [disp(dim), ..., lagrange(dim)]
 * Auxiliary fields = None
 */
void
pylith::fekernels::FaultCohesiveKin::Jf0ul_pos(const PylithInt dim,
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
                                               PylithScalar Jf0[]) {
    assert(numS >= 2);
    assert(Jf0);
    assert(sOff);
    assert(n);

    const PylithInt spaceDim = dim + 1; // :KLUDGE: dim passed in is spaceDim-1

    const PylithInt ncols = spaceDim;

    for (PylithInt i = 0; i < spaceDim; ++i) {
        Jf0[i*ncols+i] += +1.0;
    } // for
} // Jg0ul_pos


// ----------------------------------------------------------------------
// Jg0 function for integration of the slip constraint equation.
void
pylith::fekernels::FaultCohesiveKin::Jf0lu(const PylithInt dim,
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
                                           PylithScalar Jf0[]) {
    assert(numS >= 2);
    assert(Jf0);
    assert(sOff);
    assert(n);

    const PylithInt spaceDim = dim+1; // :KLUDGE: dim passed in is spaceDim-1

    const PylithInt gOffN = 0;
    const PylithInt gOffP = gOffN+spaceDim*spaceDim;
    const PylithInt ncols = spaceDim;

    for (PylithInt i = 0; i < spaceDim; ++i) {
        Jf0[gOffN+i*ncols+i] += -1.0;
        Jf0[gOffP+i*ncols+i] += +1.0;
    } // for
} // Jg0lu


// End of file
