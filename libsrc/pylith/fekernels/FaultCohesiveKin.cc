/* -*- C++ -*-
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
 * Copyright (c) 2010-2015 University of California, Davis
 *
 * See COPYING for license information.
 *
 * ----------------------------------------------------------------------
 */

#include <portinfo>

#include "pylith/fekernels/FaultCohesiveKin.hh"

#include <cassert> // USES assert()

/* ======================================================================
 * Kernels for time evolution equation with displacement and velocity
 * solution fields.
 *
 * Solution fields = [disp(dim), vel(dim), ...]
 * Auxiliary fields = None
 *
 * \int_V \vec{\phi}_v \cdot \left( \frac{\partial \vec{u}(t)}{\partial t} \right) \, dV =
 *   \int_V \vec{\phi}_v \cdot \vec{v}(t) \, dV.
 *
 * ======================================================================
 */

namespace pylith {
    namespace fekernels {
        class _FaultCohesiveKin {
public:

            /** Get offset in s where Lagrange multiplier field starts.
             *
             * Normally this would be sOff, but sOff doesn't account for having DOF for the two sides of the fault
             * passed to the hybrid kernels. This functions computes the correct offset into s for the Lagrange
             * multiplier field.
             *
             * @param[in] sOff Offset of registered subfields in solution field [numS].
             * @param[in] numS Number of registered subfields in solution field.
             *
             * @returns Offset of Lagrange multiplier field in s.
             */
            static PylithInt lagrange_sOff(const PylithInt sOff[],
                                           const PylithInt numS);

            /** Get offset in residual where Lagrange multiplier field starts.
             *
             * Normally it would be zero, but the Lagrange multiplier field is offset.
             *
             * @param[in] sOff Offset of registered subfields in solution field [numS].
             * @param[in] numS Number of registered subfields in solution field.
             *
             * @returns Offset of Lagrange multiplier field in residual.
             */
            static PylithInt lagrange_rOff(const PylithInt sOff[],
                                           const PylithInt numS);

            /* Compute tangential directions for 3-D fault.
             *
             * @param[in] dim Spatial dimension.
             * @param[in] refDir1 First choice for reference direction.
             * @param[in] refDir2 Second choice for reference direction if first fails.
             * @param[in] normDir Normal direction.
             * @param[out] tanDir1 First tangential direction.
             * @param[out] tanDIr2 Second tangential direction.
             */
            static void tangential_directions(const PylithInt dim,
                                              const PylithScalar refDir1[],
                                              const PylithScalar refDir2[],
                                              const PylithScalar normDir[],
                                              PylithScalar tanDir1[],
                                              PylithScalar tanDir2[]);

        }; // _FaultCohesiveKin
    } // fekernels
} // pylith

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
// Get offset in residual where Lagrange multiplier field starts.
PylithInt
pylith::fekernels::_FaultCohesiveKin::lagrange_rOff(const PylithInt sOff[],
                                                    const PylithInt numS) {
    PylithInt off = 0;
    const PylithInt numCount = numS - 1; // Don't include last field (Lagrange multiplier)
    for (PylithInt i = 0; i < numCount; ++i) {
        off += (sOff[i+1] - sOff[i]);
    } // for
    return off;
} // lagrange_rOff


// ----------------------------------------------------------------------
// Compute tangential directions from reference direction (first and second choice) and normal direction in 3-D.
void
pylith::fekernels::_FaultCohesiveKin::tangential_directions(const PylithInt dim,
                                                            const PylithScalar refDir1[],
                                                            const PylithScalar refDir2[],
                                                            const PylithScalar normDir[],
                                                            PylithScalar tanDir1[],
                                                            PylithScalar tanDir2[]) {
    assert(3 == dim);
    assert(refDir1);
    assert(refDir2);
    assert(normDir);
    assert(tanDir1);
    assert(tanDir2);

    const PylithInt _dim = 3;
    PylithScalar refDir[3] = { refDir1[0], refDir1[1], refDir1[2] };
    if (fabs(refDir[0]*normDir[0] + refDir[1]*normDir[1] + refDir[2]*normDir[2]) > 0.98) {
        for (PylithInt i = 0; i < _dim; ++i) {
            refDir[i] = refDir2[i];
        } // for
    } // if

    // refDir x normDir
    tanDir1[0] = +refDir[1]*normDir[2] - refDir[2]*normDir[1];
    tanDir1[1] = +refDir[2]*normDir[0] - refDir[0]*normDir[2];
    tanDir1[2] = +refDir[0]*normDir[1] - refDir[1]*normDir[0];

    // normDir x tanDir1
    tanDir2[0] = +normDir[1]*tanDir1[2] - normDir[2]*tanDir1[1];
    tanDir2[1] = +normDir[2]*tanDir1[0] - normDir[0]*tanDir1[2];
    tanDir2[2] = +normDir[0]*tanDir1[1] - normDir[1]*tanDir1[0];
} // _tangential_directions


// ----------------------------------------------------------------------
// g0 function for integration of the elasticity equation: g0u = -\lambda.
void
pylith::fekernels::FaultCohesiveKin::g0u(const PylithInt dim,
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
                                         PylithScalar g0[]) {
    assert(sOff);
    assert(s);
    assert(g0);

    assert(numS >= 2);

    const PylithInt spaceDim = dim + 1; // :KLUDGE: dim passed in is spaceDim-1

    const PylithInt gOffDispN = 0;
    const PylithInt gOffDispP = gOffDispN + spaceDim;
    const PylithInt sOffLagrange = pylith::fekernels::_FaultCohesiveKin::lagrange_sOff(sOff, numS);
    const PylithScalar* lagrange = &s[sOffLagrange];

    for (PylithInt i = 0; i < spaceDim; ++i) {
        g0[gOffDispN+i] += +lagrange[i];
        g0[gOffDispP+i] += -lagrange[i];
    } // for
} // g0u


// ----------------------------------------------------------------------
// g0 function for integration of the slip constraint equation: g0\lambda = d - u^+ + u^-.
void
pylith::fekernels::FaultCohesiveKin::g0l(const PylithInt dim,
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
                                         PylithScalar g0[]) {
    assert(sOff);
    assert(aOff);
    assert(s);
    assert(a);
    assert(g0);

    assert(numS >= 2);
    assert(numA >= 1);

    const PylithInt spaceDim = dim + 1; // :KLUDGE: dim passed in is spaceDim-1
    const PylithInt i_slip = numA-1;
    const PylithInt i_disp = 0;

    const PylithScalar* slip = &a[aOff[i_slip]];

    const PylithInt sOffDispN = sOff[i_disp];
    const PylithInt sOffDispP = sOff[i_disp]+spaceDim;
    const PylithInt gOffLagrange = 0;//pylith::fekernels::_FaultCohesiveKin::lagrange_rOff(sOff, numS);

    const PylithScalar* dispN = &s[sOffDispN];
    const PylithScalar* dispP = &s[sOffDispP];

    switch (spaceDim) {
    case 2: {
        const PylithInt _spaceDim = 2;
        const PylithScalar tanDir[2] = {-n[1], n[0] };
        for (PylithInt i = 0; i < _spaceDim; ++i) {
            const PylithScalar slipXY = n[i]*slip[0] + tanDir[i]*slip[1];
            g0[gOffLagrange+i] += slipXY - dispP[i] + dispN[i];
        } // for
        break;
    } // case 2
    case 3: {
        const PylithInt _spaceDim = 3;
        const PylithScalar* refDir1 = &constants[0];
        const PylithScalar* refDir2 = &constants[3];
        PylithScalar tanDir1[3], tanDir2[3];
        pylith::fekernels::_FaultCohesiveKin::tangential_directions(_spaceDim, refDir1, refDir2, n, tanDir1, tanDir2);

        for (PylithInt i = 0; i < _spaceDim; ++i) {
            const PylithScalar slipXYZ = n[i]*slip[0] + tanDir1[i]*slip[1] + tanDir2[i]*slip[2];
            g0[gOffLagrange+i] += slipXYZ - dispP[i] + dispN[i];
        } // for
        break;
    } // case 3
    default:
        assert(0);
    } // switch
} // g0l


// ----------------------------------------------------------------------
/* Jg0 function for integration of the displacement equation.
 *
 * Solution fields = [disp(dim), ..., lagrange(dim)]
 * Auxiliary fields = None
 */
void
pylith::fekernels::FaultCohesiveKin::Jg0ul(const PylithInt dim,
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
                                           PylithScalar Jg0[]) {
    assert(numS >= 2);
    assert(Jg0);
    assert(sOff);
    assert(n);

    const PylithInt spaceDim = dim + 1; // :KLUDGE: dim passed in is spaceDim-1

    const PylithInt gOffDispN = 0;
    const PylithInt gOffDispP = 0+spaceDim;
    const PylithInt ncols = spaceDim;

    for (PylithInt i = 0; i < spaceDim; ++i) {
        Jg0[(gOffDispN+i)*ncols+i] += +1.0;
        Jg0[(gOffDispP+i)*ncols+i] += -1.0;
    } // for
} // Jg0ul


// ----------------------------------------------------------------------
/* Jg0 function for integration of the slip constraint equation.
 *
 * Solution fields = [disp(dim), ..., lagrange(dim)]
 * Auxiliary fields = None
 */
void
pylith::fekernels::FaultCohesiveKin::Jg0lu(const PylithInt dim,
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
                                           PylithScalar Jg0[]) {
    assert(numS >= 2);
    assert(Jg0);
    assert(sOff);
    assert(n);

    const PylithInt spaceDim = dim+1; // :KLUDGE: dim passed in is spaceDim-1

    const PylithInt gOffDispN = 0;
    const PylithInt gOffDispP = 0+spaceDim;
    const PylithInt ncols = 2*spaceDim;

    for (PylithInt i = 0; i < spaceDim; ++i) {
        Jg0[i*ncols+gOffDispN+i] += +1.0;
        Jg0[i*ncols+gOffDispP+i] += -1.0;
    } // for
} // Jg0lu


// End of file
