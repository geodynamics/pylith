// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University at Buffalo
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2021 University of California, Davis
//
// See LICENSE.md for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "KinSrcLiuCos.hh" // implementation of object methods

#include "pylith/faults/KinSrcAuxiliaryFactory.hh" // USES KinSrcAuxiliaryFactory

#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_*
#include "pylith/utils/error.hh" // USES PYLITH_METHOD_BEGIN/END

#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys

#include <cassert> // USES assert()

// ---------------------------------------------------------------------------------------------------------------------
// Default constructor.
pylith::faults::KinSrcLiuCos::KinSrcLiuCos(void) {
    pylith::utils::PyreComponent::setName("kinsrcliucos");
} // constructor


// ---------------------------------------------------------------------------------------------------------------------
// Destructor.
pylith::faults::KinSrcLiuCos::~KinSrcLiuCos(void) {}


// ---------------------------------------------------------------------------------------------------------------------
// Slip time function kernel.
void
pylith::faults::KinSrcLiuCos::slipFn(const PylithInt dim,
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
                                     PylithScalar slip[]) {
    const PylithInt _numA = 3;

    assert(_numA == numA);
    assert(aOff);
    assert(a);
    assert(slip);

    const PylithInt i_initiationTime = 0;
    const PylithInt i_finalSlip = 1;
    const PylithInt i_riseTime = 2;
    const PylithScalar initiationTime = a[aOff[i_initiationTime]];
    const PylithScalar* finalSlip = &a[aOff[i_finalSlip]];
    const PylithScalar riseTime = a[aOff[i_riseTime]];

    const PylithInt i_originTime = 0;
    const PylithScalar originTime = constants[i_originTime];
    const PylithScalar t0 = originTime + initiationTime;

    if (t >= t0) {
        const PylithScalar tR = t - t0;
        const PylithScalar tau = riseTime * 1.525;
        const PylithScalar tau1 = 0.13 * tau;
        const PylithScalar tau2 = tau - tau1;
        const PylithScalar Cn = M_PI /  (1.4 * M_PI * tau1 + 1.2 * tau1 + 0.3 * M_PI * tau2);

        PylithScalar slipAmplitude = 0.0;
        if (tR <= tau1) {
            slipAmplitude = 0.7*tR - 0.7*tau1/M_PI*sin(M_PI*tR/tau1) - 0.6*tau1/(0.5*M_PI)*(cos(0.5*M_PI*tR/tau1) - 1.0);
            slipAmplitude *= Cn;
        } else if (tR <= 2.0*tau1) {
            slipAmplitude = 1.0*tR - 0.7*tau1/M_PI*sin(M_PI*tR/tau1) + 0.3*tau2/M_PI*sin(M_PI*(tR-tau1)/tau2) + 1.2*tau1/M_PI
                            - 0.3*tau1;
            slipAmplitude *= Cn;
        } else if (tR <= tau) {
            slipAmplitude = 0.3*tR + 0.3*tau2/M_PI*sin(M_PI*(tR-tau1)/tau2) + 1.1*tau1 + 1.2*tau1/M_PI;
            slipAmplitude *= Cn;
        } else {
            slipAmplitude = 1.0;
        } // if
        for (PylithInt i = 0; i < dim; ++i) {
            slip[i] = finalSlip[i] * slipAmplitude;
        } // for
    } // if
} // slipFn


// ---------------------------------------------------------------------------------------------------------------------
// Slip rate time function kernel.
void
pylith::faults::KinSrcLiuCos::slipRateFn(const PylithInt dim,
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
                                         PylithScalar slipRate[]) {
    const PylithInt _numA = 3;

    assert(_numA == numA);
    assert(aOff);
    assert(a);
    assert(slipRate);

    const PylithInt i_initiationTime = 0;
    const PylithInt i_finalSlip = 1;
    const PylithInt i_riseTime = 2;
    const PylithScalar initiationTime = a[aOff[i_initiationTime]];
    const PylithScalar* finalSlip = &a[aOff[i_finalSlip]];
    const PylithScalar riseTime = a[aOff[i_riseTime]];

    const PylithInt i_originTime = 0;
    const PylithScalar originTime = constants[i_originTime];
    const PylithScalar t0 = originTime + initiationTime;

    if (t >= t0) {
        const PylithScalar tR = t - t0;
        const PylithScalar tau = riseTime * 1.525;
        const PylithScalar tau1 = 0.13 * tau;
        const PylithScalar tau2 = tau - tau1;
        const PylithScalar Cn = M_PI /  (1.4 * M_PI * tau1 + 1.2 * tau1 + 0.3 * M_PI * tau2);

        PylithScalar slipRateAmplitude = 0.0;
        if (tR <= tau1) {
            slipRateAmplitude = 0.7 - 0.7*cos(M_PI*tR/tau1) + 0.6*sin(0.5*M_PI*tR/tau1);
            slipRateAmplitude *= Cn;
        } else if (tR <= 2.0*tau1) {
            slipRateAmplitude = 1.0 - 0.7*cos(M_PI*tR/tau1) + 0.3*cos(M_PI*(tR-tau1)/tau2);
            slipRateAmplitude *= Cn;
        } else if (tR <= tau) {
            slipRateAmplitude = 0.3 + 0.3*cos(M_PI*(tR-tau1)/tau2);
            slipRateAmplitude *= Cn;
        } else {
            slipRateAmplitude = 0.0;
        } // if
        for (PylithInt i = 0; i < dim; ++i) {
            slipRate[i] = finalSlip[i] * slipRateAmplitude;
        } // for
    } // if
} // slipRateFn


// ---------------------------------------------------------------------------------------------------------------------
// Preinitialize earthquake source. Set names/sizes of auxiliary subfields.
void
pylith::faults::KinSrcLiuCos::_auxiliaryFieldSetup(const spatialdata::units::Nondimensional& normalizer,
                                                   const spatialdata::geocoords::CoordSys* cs) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("_auxiliaryFieldSetup()");

    assert(_auxiliaryFactory);
    assert(cs);
    _auxiliaryFactory->initialize(_auxiliaryField, normalizer, cs->getSpaceDim());

    // :ATTENTION: The order for adding subfields must match the order of the auxiliary fields in the slip time function
    // kernel.

    _auxiliaryFactory->addInitiationTime(); // 0
    _auxiliaryFactory->addFinalSlip(); // 1
    _auxiliaryFactory->addRiseTime(); // 2

    _slipFnKernel = pylith::faults::KinSrcLiuCos::slipFn;
    _slipRateFnKernel = pylith::faults::KinSrcLiuCos::slipRateFn;

    PYLITH_METHOD_END;
} // _auxiliaryFieldSetup


// End of file
