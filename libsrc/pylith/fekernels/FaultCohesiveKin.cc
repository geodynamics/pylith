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

namespace pylith {
    namespace fekernels {
        class _FaultCohesiveKin {
public:

            /** Get offset in s where velocity subfield starts.
             *
             * Normally this would be sOff, but sOff doesn't account for having DOF for the two sides of the fault
             * passed to the hybrid kernels. This functions computes the correct offset into s for the velocity
             * subfield.
             *
             * @param[in] sOff Offset of registered subfields in solution field [numS].
             * @param[in] numS Number of registered subfields in solution field.
             *
             * @returns Offset of velocity subfield in s.
             */
            static PylithInt velocity_sOff(const PylithInt sOff[],
                                           const PylithInt numS);

            /** Get offset in s where Lagrange multiplier subfield starts.
             *
             * Normally this would be sOff, but sOff doesn't account for having DOF for the two sides of the fault
             * passed to the hybrid kernels. This functions computes the correct offset into s for the Lagrange
             * multiplier subfield.
             *
             * @param[in] sOff Offset of registered subfields in solution field [numS].
             * @param[in] numS Number of registered subfields in solution field.
             *
             * @returns Offset of Lagrange multiplier subfield in s.
             */
            static PylithInt lagrange_sOff(const PylithInt sOff[],
                                           const PylithInt numS);

        }; // _FaultCohesiveKin
    } // fekernels
} // pylith

// ----------------------------------------------------------------------
// Get offset in s where velocity subfield starts.
PylithInt
pylith::fekernels::_FaultCohesiveKin::velocity_sOff(const PylithInt sOff[],
                                                    const PylithInt numS) {
    PylithInt off = 0;
    const PylithInt numCount = 1; // [displacement, velocity, ...]
    for (PylithInt i = 0; i < numCount; ++i) {
        off += 2*(sOff[i+1] - sOff[i]);
    } // for
    return off;
} // velocity_sOff


// ----------------------------------------------------------------------
// Get offset in s where Lagrange multiplier field starts.
PylithInt
pylith::fekernels::_FaultCohesiveKin::lagrange_sOff(const PylithInt sOff[],
                                                    const PylithInt numS) {
    PylithInt off = 0;
    const PylithInt numCount = numS - 1; // Don't include last field (Lagrange multiplier)
    for (PylithInt i = 0; i < numCount; ++i) {
        off += 2*(sOff[i+1] - sOff[i]);
    } // for
    return off;
} // lagrange_sOff


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

    const PylithInt fOffN = 0;
    const PylithInt sOffLagrange = pylith::fekernels::_FaultCohesiveKin::lagrange_sOff(sOff, numS);
    const PylithScalar* lagrange = &s[sOffLagrange];

    for (PylithInt i = 0; i < spaceDim; ++i) {
        f0[fOffN+i] += -lagrange[i];
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

    const PylithInt fOffP = 0;
    const PylithInt sOffLagrange = pylith::fekernels::_FaultCohesiveKin::lagrange_sOff(sOff, numS);
    const PylithScalar* lagrange = &s[sOffLagrange];

    for (PylithInt i = 0; i < spaceDim; ++i) {
        // f0[fOffN+i] += -lagrange[i];
        f0[fOffP+i] += +lagrange[i];
    } // for
} // f0u_pos


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
    const PylithInt i_slip = numA-1;
    const PylithInt i_disp = 0;

    const PylithScalar* slip = &a[aOff[i_slip]];

    const PylithInt sOffDispN = sOff[i_disp];
    const PylithInt sOffDispP = sOffDispN+spaceDim;
    const PylithInt fOffLagrange = 0;

    const PylithScalar* dispN = &s[sOffDispN];
    const PylithScalar* dispP = &s[sOffDispP];

    switch (spaceDim) {
    case 2: {
        const PylithInt _spaceDim = 2;
        const PylithScalar tanDir[2] = {-n[1], n[0] };
        for (PylithInt i = 0; i < _spaceDim; ++i) {
            const PylithScalar slipXY = n[i]*slip[0] + tanDir[i]*slip[1];
            f0[fOffLagrange+i] += dispP[i] - dispN[i] - slipXY;
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
            const PylithScalar slipXYZ = n[i]*slip[0] + tanDir1[i]*slip[1] + tanDir2[i]*slip[2];
            f0[fOffLagrange+i] += dispP[i] - dispN[i] - slipXYZ;
        } // for
        break;
    } // case 3
    default:
        assert(0);
    } // switch
} // f0l_u


// ----------------------------------------------------------------------
// f0 function for slip rate constraint equation: f0\lambda = (v^+ - v^-) - \dot{d}
void
pylith::fekernels::FaultCohesiveKin::f0l_v(const PylithInt dim,
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
    const PylithInt i_slipRate = numA-1;

    const PylithScalar* slipRate = &a[aOff[i_slipRate]];

    const PylithInt sOffVelN = _FaultCohesiveKin::velocity_sOff(sOff, numS);
    const PylithInt sOffVelP = sOffVelN+spaceDim;
    const PylithInt fOffLagrange = 0;

    const PylithScalar* velN = &s[sOffVelN];
    const PylithScalar* velP = &s[sOffVelP];

    switch (spaceDim) {
    case 2: {
        const PylithInt _spaceDim = 2;
        const PylithScalar tanDir[2] = {-n[1], n[0] };
        for (PylithInt i = 0; i < _spaceDim; ++i) {
            const PylithScalar slipRateXY = n[i]*slipRate[0] + tanDir[i]*slipRate[1];
            f0[fOffLagrange+i] += velP[i] - velN[i] - slipRateXY;
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
            const PylithScalar slipRateXYZ = n[i]*slipRate[0] + tanDir1[i]*slipRate[1] + tanDir2[i]*slipRate[2];
            f0[fOffLagrange+i] += velP[i] - velN[i] - slipRateXYZ;
        } // for
        break;
    } // case 3
    default:
        assert(0);
    } // switch
} // f0l_v


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

    const PylithInt sOffVelN = _FaultCohesiveKin::velocity_sOff(sOff, numS);
    const PylithInt sOffVelP = sOffVelN+spaceDim;
    const PylithInt fOffLagrange = 0;

    const PylithScalar* accN = &s_t[sOffVelN];
    const PylithScalar* accP = &s_t[sOffVelP];

    switch (spaceDim) {
    case 2: {
        const PylithInt _spaceDim = 2;
        const PylithScalar tanDir[2] = {-n[1], n[0] };
        for (PylithInt i = 0; i < _spaceDim; ++i) {
            const PylithScalar slipAccXY = n[i]*slipAcc[0] + tanDir[i]*slipAcc[1];
            f0[fOffLagrange+i] += accP[i] - accN[i] - slipAccXY;
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
            f0[fOffLagrange+i] += accP[i] - accN[i] - slipAccXYZ;
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

    const PylithInt gOffN = 0;
    const PylithInt ncols = spaceDim;

    for (PylithInt i = 0; i < spaceDim; ++i) {
        Jf0[(gOffN+i)*ncols+i] += -1.0;
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

    const PylithInt gOffN = 0;
    const PylithInt gOffP = gOffN+spaceDim;
    const PylithInt ncols = spaceDim;

    for (PylithInt i = 0; i < spaceDim; ++i) {
        Jf0[(gOffP+i)*ncols+i] += +1.0;
    } // for
} // Jg0ul_pos


// ----------------------------------------------------------------------
/* Jg0 function for integration of the slip constraint equation.
 *
 * Solution fields = [disp(dim), ..., lagrange(dim)]
 * Auxiliary fields = None
 */
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
    const PylithInt gOffP = gOffN+spaceDim;
    const PylithInt ncols = 2*spaceDim;

    for (PylithInt i = 0; i < spaceDim; ++i) {
        Jf0[i*ncols+gOffN+i] += -1.0;
        Jf0[i*ncols+gOffP+i] += +1.0;
    } // for
} // Jg0lu


// End of file
