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
    const PylithInt gOffLagrange = pylith::fekernels::_FaultCohesiveKin::lagrange_rOff(sOff, numS);

    const PylithScalar* dispN = &s[sOffDispN];
    const PylithScalar* dispP = &s[sOffDispP];

    for (PylithInt i = 0; i < spaceDim; ++i) {
        g0[gOffLagrange+i] += slip[i] - dispP[i] + dispN[i];
    } // for
} // g0l


#include <iostream>
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

    const PylithInt i_disp = 0;
    const PylithInt sOffDispN = sOff[i_disp];
    const PylithInt sOffDispP = sOff[i_disp]+spaceDim;
    const PylithInt sOffLagrange = 0;
    const PylithInt ncols = sOffLagrange + spaceDim;

    for (PylithInt i = 0; i < spaceDim; ++i) {
        std::cout << "Jg0ul ncols="<<ncols<<", uN: r="<<(sOffDispN+i)<<", c="<<sOffLagrange+i<<", index="<<(sOffDispN+i)*ncols+sOffLagrange+i<<std::endl;
        std::cout << "Jg0ul ncols="<<ncols<<", uP: r="<<(sOffDispP+i)<<", c="<<sOffLagrange+i<<", index="<<(sOffDispP+i)*ncols+sOffLagrange+i<<std::endl;
        Jg0[(sOffDispN+i)*ncols+sOffLagrange+i] += +1.0;
        Jg0[(sOffDispP+i)*ncols+sOffLagrange+i] += -1.0;
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

    const PylithInt i_disp = 0;
    const PylithInt sOffDispN = sOff[i_disp];
    const PylithInt sOffDispP = sOff[i_disp]+spaceDim;
    const PylithInt sOffLagrange = 0; // pylith::fekernels::_FaultCohesiveKin::lagrange_sOff(sOff, numS);
    const PylithInt ncols = 2*spaceDim;

    for (PylithInt i = 0; i < spaceDim; ++i) {
        std::cout << "Jg0lu uN: r="<<0*sOffLagrange+i<<", c="<<(sOffDispN+i)<<", index="<<(0*sOffLagrange+i)*ncols+sOffDispN+i<<std::endl;
        std::cout << "Jg0lu uP: r="<<0*sOffLagrange+i<<", c="<<(sOffDispP+i)<<", index="<<(0*sOffLagrange+i)*ncols+sOffDispP+i<<std::endl;
        Jg0[(0*sOffLagrange+i)*ncols+sOffDispN+i] += +1.0;
        Jg0[(0*sOffLagrange+i)*ncols+sOffDispP+i] += -1.0;
    } // for
} // Jg0lu


// End of file
