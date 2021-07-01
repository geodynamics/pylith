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

#include <portinfo>

#include "pylith/fekernels/NeumannTimeDependent.hh"

#include <cassert> // USES assert()

/* ======================================================================
 * Kernels for computing value from parameters for time-dependent boundary conditions.
 * ======================================================================
 */

/* ----------------------------------------------------------------------
 * Scalar initial value term for time-dependent boundary condition.
 *
 * f_0(x)
 */
void
pylith::fekernels::NeumannTimeDependent::g0_initial_scalar(const PylithInt dim,
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
    const PylithInt _numA = 1;
    assert(_numA == numA);
    assert(aOff);
    assert(a);
    assert(g0);

    const PylithInt i_initial = aOff[0];

    g0[0] += a[i_initial];
} // g0_initial_scalar


/* ----------------------------------------------------------------------
 * Vector initial value term for time-dependent boundary condition.
 *
 * f_0(x)
 */
void
pylith::fekernels::NeumannTimeDependent::g0_initial_vector(const PylithInt dim,
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
    const PylithInt _numA = 1;
    assert(_numA == numA);
    assert(aOff);
    assert(a);
    assert(6 == numConstants);
    assert(constants);
    assert(2 == dim || 3 == dim);
    assert(g0);

    const PylithInt i_initial = aOff[0];

    switch (dim) {
    case 2: {
        const PylithInt _dim = 2;
        const PylithScalar tanDir[2] = {-n[1], n[0] };
        for (PylithInt i = 0; i < _dim; ++i) {
            g0[i] += a[i_initial+0]*tanDir[i] + a[i_initial+1]*n[i];
        } // for
        break;
    } // case 2
    case 3: {
        const PylithInt _dim = 3;
        const PylithScalar* refDir1 = &constants[0];
        const PylithScalar* refDir2 = &constants[3];
        PylithScalar tanDir1[3], tanDir2[3];
        _tangential_directions(_dim, refDir1, refDir2, n, tanDir1, tanDir2);

        for (PylithInt i = 0; i < _dim; ++i) {
            g0[i] += a[i_initial+0]*tanDir1[i] + a[i_initial+1]*tanDir2[i] + a[i_initial+2]*n[i];
        } // for
        break;
    } // case 3
    default:
        assert(0);
    } // switch

} // g0_initial_vector


/* ----------------------------------------------------------------------
 * Scalar rate term for time-dependent boundary condition.
 *
 * \dot{f}_1(x) * (t-t_1(x)) for t >= t_1(x).
 */
void
pylith::fekernels::NeumannTimeDependent::g0_rate_scalar(const PylithInt dim,
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
    const PylithInt _numA = 2;
    assert(_numA == numA);
    assert(aOff);
    assert(a);
    assert(g0);

    const PylithInt i_rate = aOff[0];
    const PylithInt i_start = aOff[1];

    const PylithScalar tRel = t - a[i_start];
    if (tRel > 0.0) {
        g0[0] += a[i_rate]*tRel;
    } // if/else

} // g0_rate_scalar


/* ----------------------------------------------------------------------
 * Vector rate term for time-dependent boundary condition.
 *
 * \dot{f}_1(x) * (t-t_1(x)) for t >= t_1(x).
 */
void
pylith::fekernels::NeumannTimeDependent::g0_rate_vector(const PylithInt dim,
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
    const PylithInt _numA = 2;
    assert(_numA == numA);
    assert(aOff);
    assert(a);
    assert(6 == numConstants);
    assert(constants);
    assert(2 == dim || 3 == dim);
    assert(g0);

    const PylithInt i_rate = aOff[0];
    const PylithInt i_start = aOff[1];

    const PylithScalar tRel = t - a[i_start];
    if (tRel > 0.0) {
        switch (dim) {
        case 2: {
            const PylithInt _dim = 2;
            const PylithScalar tanDir[2] = {-n[1], n[0] };
            for (PylithInt i = 0; i < _dim; ++i) {
                g0[i] += tRel * (a[i_rate+0]*tanDir[i] + a[i_rate+1]*n[i]);
            } // for
            break;
        } // case 2
        case 3: {
            const PylithInt _dim = 3;
            const PylithScalar* refDir1 = &constants[0];
            const PylithScalar* refDir2 = &constants[3];
            PylithScalar tanDir1[3], tanDir2[3];
            _tangential_directions(_dim, refDir1, refDir2, n, tanDir1, tanDir2);

            for (PylithInt i = 0; i < _dim; ++i) {
                g0[i] += tRel * (a[i_rate+0]*tanDir1[i] + a[i_rate+1]*tanDir2[i] + a[i_rate+2]*n[i]);
            } // for
            break;
        } // case 3
        } // switch
    } // if
} // g0_rate_vector


/* ----------------------------------------------------------------------
 * Scalar time history term for time-dependent boundary condition.
 *
 * f_2(x) * a(t-t_2(x)) for t >= t_2(x).
 */
void
pylith::fekernels::NeumannTimeDependent::g0_timeHistory_scalar(const PylithInt dim,
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
    const PylithInt _numA = 3;
    assert(_numA == numA);
    assert(aOff);
    assert(a);
    assert(g0);

    const PylithInt i_amplitude = aOff[0];
    const PylithInt i_value = aOff[2];

    g0[0] += a[i_amplitude]*a[i_value];
} // g0_timeHistory_scalar


/* ----------------------------------------------------------------------
 * Vector time history term for time-dependent boundary condition.
 *
 * f_2(x) * a(t-t_2(x)) for t >= t_2(x).
 */
void
pylith::fekernels::NeumannTimeDependent::g0_timeHistory_vector(const PylithInt dim,
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
    const PylithInt _numA = 3;
    assert(_numA == numA);
    assert(aOff);
    assert(a);
    assert(6 == numConstants);
    assert(constants);
    assert(2 == dim || 3 == dim);
    assert(g0);

    const PylithInt i_amplitude = aOff[0];
    const PylithInt i_value = aOff[2];

    switch (dim) {
    case 2: {
        const PylithInt _dim = 2;
        const PylithScalar tanDir[2] = {-n[1], n[0] };
        for (PylithInt i = 0; i < _dim; ++i) {
            g0[i] += a[i_value] * (a[i_amplitude+0]*tanDir[i] + a[i_amplitude+1]*n[i]);
        } // for
        break;
    } // case 2
    case 3: {
        const PylithInt _dim = 3;
        const PylithScalar* refDir1 = &constants[0];
        const PylithScalar* refDir2 = &constants[3];
        PylithScalar tanDir1[3], tanDir2[3];
        _tangential_directions(_dim, refDir1, refDir2, n, tanDir1, tanDir2);

        for (PylithInt i = 0; i < _dim; ++i) {
            g0[i] += a[i_value] * (a[i_amplitude+0]*tanDir1[i] + a[i_amplitude+1]*tanDir2[i] + a[i_amplitude+2]*n[i]);
        } // for
        break;
    } // case 3
    } // switch
} // g0_timeHistory_vector


// ----------------------------------------------------------------------
// Compute boundary condition scalar value using initial and rate terms.
void
pylith::fekernels::NeumannTimeDependent::g0_initialRate_scalar(const PylithInt dim,
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
    const PylithInt _numA = 3;
    assert(_numA == numA);
    assert(aOff);
    assert(a);
    assert(g0);

    const PylithInt i_initial = 0;
    const PylithInt i_rate = 1;
    const PylithInt i_start = 2;

    const PylithInt numAInitial = 1;
    const PylithInt aOffInitial[1] = { aOff[i_initial] };
    const PylithInt* aOffInitial_x = NULL;

    const PylithInt numARate = 2;
    const PylithInt aOffRate[2] = { aOff[i_rate], aOff[i_start] };
    const PylithInt* aOffRate_x = NULL;

    g0_initial_scalar(dim, numS, numAInitial, sOff, sOff_x, s, s_t, s_x, aOffInitial, aOffInitial_x, a, a_t, a_x,
                      t, x, n, numConstants, constants, g0);
    g0_rate_scalar(dim, numS, numARate, sOff, sOff_x, s, s_t, s_x, aOffRate, aOffRate_x, a, a_t, a_x,
                   t, x, n, numConstants, constants, g0);
} // g0_initialRate_scalar


// ----------------------------------------------------------------------
// Compute boundary condition vector value using initial and rate terms.
void
pylith::fekernels::NeumannTimeDependent::g0_initialRate_vector(const PylithInt dim,
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
    const PylithInt _numA = 3;
    assert(_numA == numA);
    assert(aOff);
    assert(a);
    assert(g0);

    const PylithInt i_initial = 0;
    const PylithInt i_rate = 1;
    const PylithInt i_start = 2;

    const PylithInt numAInitial = 1;
    const PylithInt aOffInitial[1] = { aOff[i_initial] };
    const PylithInt* aOffInitial_x = NULL;

    const PylithInt numARate = 2;
    const PylithInt aOffRate[2] = { aOff[i_rate], aOff[i_start] };
    const PylithInt* aOffRate_x = NULL;

    g0_initial_vector(dim, numS, numAInitial, sOff, sOff_x, s, s_t, s_x, aOffInitial, aOffInitial_x, a, a_t, a_x,
                      t, x, n, numConstants, constants, g0);
    g0_rate_vector(dim, numS, numARate, sOff, sOff_x, s, s_t, s_x, aOffRate, aOffRate_x, a, a_t, a_x,
                   t, x, n, numConstants, constants, g0);
} // g0_initialRate_vector


// ----------------------------------------------------------------------
// Compute boundary condition scalar value using initial and time history terms.
void
pylith::fekernels::NeumannTimeDependent::g0_initialTimeHistory_scalar(const PylithInt dim,
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
    const PylithInt _numA = 4;
    assert(_numA == numA);
    assert(aOff);
    assert(a);
    assert(g0);

    const PylithInt i_initial = 0;
    const PylithInt i_thAmp = 1;
    const PylithInt i_thStart = 2;
    const PylithInt i_thValue = 3;

    const PylithInt numAInitial = 1;
    const PylithInt aOffInitial[1] = { aOff[i_initial] };
    const PylithInt* aOffInitial_x = NULL;

    const PylithInt numATimeHistory = 3;
    const PylithInt aOffTimeHistory[3] = { aOff[i_thAmp], aOff[i_thStart], aOff[i_thValue] };
    const PylithInt* aOffTimeHistory_x = NULL;

    g0_initial_scalar(dim, numS, numAInitial, sOff, sOff_x, s, s_t, s_x, aOffInitial, aOffInitial_x, a, a_t, a_x,
                      t, x, n, numConstants, constants, g0);
    g0_timeHistory_scalar(dim, numS, numATimeHistory, sOff, sOff_x, s, s_t, s_x,
                          aOffTimeHistory, aOffTimeHistory_x, a, a_t, a_x, t, x, n, numConstants, constants, g0);
} // g0_initialTimeHistory_scalar


// ----------------------------------------------------------------------
// Compute boundary condition vector value using initial and time history terms.
void
pylith::fekernels::NeumannTimeDependent::g0_initialTimeHistory_vector(const PylithInt dim,
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
    const PylithInt _numA = 4;
    assert(_numA == numA);
    assert(aOff);
    assert(a);
    assert(g0);

    const PylithInt i_initial = 0;
    const PylithInt i_thAmp = 1;
    const PylithInt i_thStart = 2;
    const PylithInt i_thValue = 3;

    const PylithInt numAInitial = 1;
    const PylithInt aOffInitial[1] = { aOff[i_initial] };
    const PylithInt* aOffInitial_x = NULL;

    const PylithInt numATimeHistory = 3;
    const PylithInt aOffTimeHistory[3] = { aOff[i_thAmp], aOff[i_thStart], aOff[i_thValue] };
    const PylithInt* aOffTimeHistory_x = NULL;

    g0_initial_vector(dim, numS, numAInitial, sOff, sOff_x, s, s_t, s_x, aOffInitial, aOffInitial_x, a, a_t, a_x,
                      t, x, n, numConstants, constants, g0);
    g0_timeHistory_vector(dim, numS, numATimeHistory, sOff, sOff_x, s, s_t, s_x,
                          aOffTimeHistory, aOffTimeHistory_x, a, a_t, a_x, t, x, n, numConstants, constants, g0);
} // g0_initialTimeHistory_vector


// ----------------------------------------------------------------------
// Compute boundary condition scalar value using rate and time history terms.
void
pylith::fekernels::NeumannTimeDependent::g0_rateTimeHistory_scalar(const PylithInt dim,
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
    const PylithInt _numA = 5;
    assert(_numA == numA);
    assert(aOff);
    assert(a);
    assert(g0);

    const PylithInt i_rateAmp = 0;
    const PylithInt i_rateStart = 0;
    const PylithInt i_thAmp = 1;
    const PylithInt i_thStart = 2;
    const PylithInt i_thValue = 3;

    const PylithInt numARate = 1;
    const PylithInt aOffRate[2] = { aOff[i_rateAmp], aOff[i_rateStart] };
    const PylithInt* aOffRate_x = NULL;

    const PylithInt numATimeHistory = 3;
    const PylithInt aOffTimeHistory[3] = { aOff[i_thAmp], aOff[i_thStart], aOff[i_thValue] };
    const PylithInt* aOffTimeHistory_x = NULL;

    g0_rate_scalar(dim, numS, numARate, sOff, sOff_x, s, s_t, s_x, aOffRate, aOffRate_x, a, a_t, a_x,
                   t, x, n, numConstants, constants, g0);
    g0_timeHistory_scalar(dim, numS, numATimeHistory, sOff, sOff_x, s, s_t, s_x,
                          aOffTimeHistory, aOffTimeHistory_x, a, a_t, a_x, t, x, n, numConstants, constants, g0);
} // g0_rateTimeHistory_scalar


// ----------------------------------------------------------------------
// Compute boundary condition vector value using rate and time history terms.
void
pylith::fekernels::NeumannTimeDependent::g0_rateTimeHistory_vector(const PylithInt dim,
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
    const PylithInt _numA = 5;
    assert(_numA == numA);
    assert(aOff);
    assert(a);
    assert(g0);

    const PylithInt i_rateAmp = 0;
    const PylithInt i_rateStart = 0;
    const PylithInt i_thAmp = 1;
    const PylithInt i_thStart = 2;
    const PylithInt i_thValue = 3;

    const PylithInt numARate = 1;
    const PylithInt aOffRate[2] = { aOff[i_rateAmp], aOff[i_rateStart] };
    const PylithInt* aOffRate_x = NULL;

    const PylithInt numATimeHistory = 3;
    const PylithInt aOffTimeHistory[3] = { aOff[i_thAmp], aOff[i_thStart], aOff[i_thValue] };
    const PylithInt* aOffTimeHistory_x = NULL;

    g0_rate_vector(dim, numS, numARate, sOff, sOff_x, s, s_t, s_x, aOffRate, aOffRate_x, a, a_t, a_x,
                   t, x, n, numConstants, constants, g0);
    g0_timeHistory_vector(dim, numS, numATimeHistory, sOff, sOff_x, s, s_t, s_x,
                          aOffTimeHistory, aOffTimeHistory_x, a, a_t, a_x, t, x, n, numConstants, constants, g0);
} // g0_rateTimeHistory_vector


// ----------------------------------------------------------------------
// Compute boundary condition scalar value using initial, rate ,and time history terms.
void
pylith::fekernels::NeumannTimeDependent::g0_initialRateTimeHistory_scalar(const PylithInt dim,
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
    const PylithInt _numA = 6;
    assert(_numA == numA);
    assert(aOff);
    assert(a);
    assert(g0);

    const PylithInt i_valueAmp = 0;
    const PylithInt i_rateAmp = 1;
    const PylithInt i_rateStart = 2;
    const PylithInt i_thAmp = 3;
    const PylithInt i_thStart = 4;
    const PylithInt i_thValue = 5;

    const PylithInt numAInitial = 1;
    const PylithInt aOffInitial[1] = { aOff[i_valueAmp] };
    const PylithInt* aOffInitial_x = NULL;

    const PylithInt numARate = 1;
    const PylithInt aOffRate[2] = { aOff[i_rateAmp], aOff[i_rateStart] };
    const PylithInt* aOffRate_x = NULL;

    const PylithInt numATimeHistory = 3;
    const PylithInt aOffTimeHistory[3] = { aOff[i_thAmp], aOff[i_thStart], aOff[i_thValue] };
    const PylithInt* aOffTimeHistory_x = NULL;

    g0_initial_scalar(dim, numS, numAInitial, sOff, sOff_x, s, s_t, s_x, aOffInitial, aOffInitial_x, a, a_t, a_x,
                      t, x, n, numConstants, constants, g0);
    g0_rate_scalar(dim, numS, numARate, sOff, sOff_x, s, s_t, s_x, aOffRate, aOffRate_x, a, a_t, a_x,
                   t, x, n, numConstants, constants, g0);
    g0_timeHistory_scalar(dim, numS, numATimeHistory, sOff, sOff_x, s, s_t, s_x,
                          aOffTimeHistory, aOffTimeHistory_x, a, a_t, a_x, t, x, n, numConstants, constants, g0);
} // g0_initialRateTimeHistory_scalar


// ----------------------------------------------------------------------
// Compute boundary condition vector value using initial, rate, and time history terms.
void
pylith::fekernels::NeumannTimeDependent::g0_initialRateTimeHistory_vector(const PylithInt dim,
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
    const PylithInt _numA = 6;
    assert(_numA == numA);
    assert(aOff);
    assert(a);
    assert(g0);

    const PylithInt i_valueAmp = 0;
    const PylithInt i_rateAmp = 1;
    const PylithInt i_rateStart = 2;
    const PylithInt i_thAmp = 3;
    const PylithInt i_thStart = 4;
    const PylithInt i_thValue = 5;

    const PylithInt numAInitial = 1;
    const PylithInt aOffInitial[1] = { aOff[i_valueAmp] };
    const PylithInt* aOffInitial_x = NULL;

    const PylithInt numARate = 1;
    const PylithInt aOffRate[2] = { aOff[i_rateAmp], aOff[i_rateStart] };
    const PylithInt* aOffRate_x = NULL;

    const PylithInt numATimeHistory = 3;
    const PylithInt aOffTimeHistory[3] = { aOff[i_thAmp], aOff[i_thStart], aOff[i_thValue] };
    const PylithInt* aOffTimeHistory_x = NULL;

    g0_initial_vector(dim, numS, numAInitial, sOff, sOff_x, s, s_t, s_x, aOffInitial, aOffInitial_x, a, a_t, a_x,
                      t, x, n, numConstants, constants, g0);
    g0_rate_vector(dim, numS, numARate, sOff, sOff_x, s, s_t, s_x, aOffRate, aOffRate_x, a, a_t, a_x,
                   t, x, n, numConstants, constants, g0);
    g0_timeHistory_vector(dim, numS, numATimeHistory, sOff, sOff_x, s, s_t, s_x,
                          aOffTimeHistory, aOffTimeHistory_x, a, a_t, a_x, t, x, n, numConstants, constants, g0);
} // g0_initialRateTimeHistory_vector


// ----------------------------------------------------------------------
// Compute tangential directions from reference direction (first and second choice) and normal direction in 3-D.
void
pylith::fekernels::NeumannTimeDependent::_tangential_directions(const PylithInt dim,
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


// End of file
