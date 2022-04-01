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

#include "KinSrcRamp.hh" // implementation of object methods

#include "pylith/faults/KinSrcAuxiliaryFactory.hh" // USES KinSrcAuxiliaryFactory

#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_*
#include "pylith/utils/error.hh" // USES PYLITH_METHOD_BEGIN/END

#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys

#include <cassert> // USES assert()

namespace pylith {
    namespace faults {
        namespace _KinSrcRamp {
            double maxAcc(const double finalSlip,
                          const double riseTime,
                          const double impulseDuration);

        } // _KinSrcRamp
    } // faults
} // pylith

// ------------------------------------------------------------------------------------------------
// Default constructor.
pylith::faults::KinSrcRamp::KinSrcRamp(void) {
    pylith::utils::PyreComponent::setName("kinsrcramp");
} // constructor


// ------------------------------------------------------------------------------------------------
// Destructor.
pylith::faults::KinSrcRamp::~KinSrcRamp(void) {}


// ------------------------------------------------------------------------------------------------
// Slip time function kernel.
void
pylith::faults::KinSrcRamp::slipFn(const PylithInt dim,
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
    const PylithInt _numA = 4;

    assert(_numA == numA);
    assert(aOff);
    assert(a);
    assert(slip);

    const PylithInt i_initiationTime = 0;
    const PylithInt i_finalSlip = 1;
    const PylithInt i_riseTime = 2;
    const PylithInt i_impulseDuration = 3;
    const PylithScalar initiationTime = a[aOff[i_initiationTime]];
    const PylithScalar* finalSlip = &a[aOff[i_finalSlip]];
    const PylithScalar riseTime = a[aOff[i_riseTime]];
    const PylithScalar impulseDuration = a[aOff[i_impulseDuration]];

    const PylithInt i_originTime = 0;
    const PylithScalar originTime = constants[i_originTime];
    const PylithScalar t0 = originTime + initiationTime;

    const double finalSlipMag = dim == 2 ?
                                sqrt(pow(finalSlip[0],2) + pow(finalSlip[1],2)) :
                                sqrt(pow(finalSlip[0],2) + pow(finalSlip[1],2) + pow(finalSlip[2],2));
    const double slipCoef = 2.0 * _KinSrcRamp::maxAcc(finalSlipMag, riseTime, impulseDuration) / impulseDuration;

    double slipTimeFn = 0.0;
    if (t-t0 < 0.0) {
        slipTimeFn = 0.0;
    } else if (t-t0 <= 0.5*impulseDuration) {
        slipTimeFn = 1/6.0*pow(t,3);
    } else if (t-t0 <= impulseDuration) {
        slipTimeFn =
            0.5*impulseDuration*pow(t,2)
            -1/6.0*pow(t,3)
            -0.25*pow(impulseDuration,2)*t
            +1.0/24.0*pow(impulseDuration,3);
    } else if (t-t0 <= riseTime-impulseDuration) {
        slipTimeFn =
            0.25*pow(impulseDuration,2)*t
            -0.125*pow(impulseDuration,3);
    } else if (t-t0 <= riseTime-0.5*impulseDuration) {
        slipTimeFn =
            -1/6.0*pow(t,3)
            +0.5*(riseTime-impulseDuration)*pow(t,2)
            -0.5*pow(riseTime-impulseDuration,2)*t
            +0.25*pow(impulseDuration,2)*t
            +1/6.0*pow(riseTime-impulseDuration,3)
            -0.125*pow(impulseDuration,3);
    } else if (t-t0 <= riseTime) {
        const double ti = riseTime - 0.5*impulseDuration;
        slipTimeFn =
            1/6.0*pow(t,3)
            -0.5*riseTime*pow(t,2)
            +0.5*pow(riseTime,2)*t
            -1/3.0*pow(ti,3)
            +0.5*riseTime*pow(ti,2)
            -0.5*pow(riseTime,2)*ti
            +0.5*(riseTime-impulseDuration)*pow(ti,2)
            -0.5*pow(riseTime-impulseDuration,2)*ti
            +0.25*pow(impulseDuration,2)*ti
            +1/6.0*pow(riseTime-impulseDuration,3)
            -0.125*pow(impulseDuration,3);
    } else {
        const double ti = riseTime - 0.5*impulseDuration;
        slipTimeFn =
            1/6.0*pow(riseTime,3)
            -1/3.0*pow(ti,3)
            +0.5*riseTime*pow(ti,2)
            -0.5*pow(riseTime,2)*ti
            +0.5*(riseTime-impulseDuration)*pow(ti,2)
            -0.5*pow(riseTime-impulseDuration,2)*ti
            +0.25*pow(impulseDuration,2)*ti
            +1/6.0*pow(riseTime-impulseDuration,3)
            -0.125*pow(impulseDuration,3);
    } // if/else

    for (PylithInt i = 0; i < dim; ++i) {
        slip[i] = slipCoef * slipTimeFn * finalSlip[i] / finalSlipMag;
    } // for
} // slipFn


// ------------------------------------------------------------------------------------------------
// Slip rate time function kernel.
void
pylith::faults::KinSrcRamp::slipRateFn(const PylithInt dim,
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
    const PylithInt _numA = 4;

    assert(_numA == numA);
    assert(aOff);
    assert(a);
    assert(slipRate);

    const PylithInt i_initiationTime = 0;
    const PylithInt i_finalSlip = 1;
    const PylithInt i_riseTime = 2;
    const PylithInt i_impulseDuration = 3;
    const PylithScalar initiationTime = a[aOff[i_initiationTime]];
    const PylithScalar* finalSlip = &a[aOff[i_finalSlip]];
    const PylithScalar riseTime = a[aOff[i_riseTime]];
    const PylithScalar impulseDuration = a[aOff[i_impulseDuration]];

    const PylithInt i_originTime = 0;
    const PylithScalar originTime = constants[i_originTime];
    const PylithScalar t0 = originTime + initiationTime;

    const double finalSlipMag = dim == 2 ?
                                sqrt(pow(finalSlip[0],2) + pow(finalSlip[1],2)) :
                                sqrt(pow(finalSlip[0],2) + pow(finalSlip[1],2) + pow(finalSlip[2],2));
    const double slipCoef = 2.0 * _KinSrcRamp::maxAcc(finalSlipMag, riseTime, impulseDuration) / impulseDuration;

    double slipRateTimeFn = 0.0;
    if (t-t0 < 0.0) {
        slipRateTimeFn = 0.0;
    } else if (t-t0 <= 0.5*impulseDuration) {
        slipRateTimeFn = 0.5*pow(t,2);
    } else if (t-t0 <= impulseDuration) {
        slipRateTimeFn = impulseDuration*t -0.5*pow(t,2) -0.25*pow(impulseDuration,2);
    } else if (t-t0 <= riseTime-impulseDuration) {
        slipRateTimeFn = 0.25*pow(impulseDuration,2);
    } else if (t-t0 <= riseTime-0.5*impulseDuration) {
        slipRateTimeFn =
            -0.5*pow(t,2)
            +(riseTime-impulseDuration)*t
            -0.5*pow(riseTime-impulseDuration,2)
            +0.25*pow(impulseDuration,2);
    } else if (t-t0 <= riseTime) {
        slipRateTimeFn = 0.5*pow(t,2) -riseTime*t +0.5*pow(riseTime,2);
    } else {
        slipRateTimeFn = 0.0;
    } // if/else

    for (PylithInt i = 0; i < dim; ++i) {
        slipRate[i] = slipCoef * slipRateTimeFn * finalSlip[i] / finalSlipMag;
    } // for
} // slipRateFn


// ------------------------------------------------------------------------------------------------
// Slip acceleration time function kernel.
void
pylith::faults::KinSrcRamp::slipAccFn(const PylithInt dim,
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
                                      PylithScalar slipAcc[]) {
    const PylithInt _numA = 4;

    assert(_numA == numA);
    assert(aOff);
    assert(a);
    assert(slipAcc);

    const PylithInt i_initiationTime = 0;
    const PylithInt i_finalSlip = 1;
    const PylithInt i_riseTime = 2;
    const PylithInt i_impulseDuration = 3;
    const PylithScalar initiationTime = a[aOff[i_initiationTime]];
    const PylithScalar* finalSlip = &a[aOff[i_finalSlip]];
    const PylithScalar riseTime = a[aOff[i_riseTime]];
    const PylithScalar impulseDuration = a[aOff[i_impulseDuration]];

    const PylithInt i_originTime = 0;
    const PylithScalar originTime = constants[i_originTime];
    const PylithScalar t0 = originTime + initiationTime;

    const double finalSlipMag = dim == 2 ?
                                sqrt(pow(finalSlip[0],2) + pow(finalSlip[1],2)) :
                                sqrt(pow(finalSlip[0],2) + pow(finalSlip[1],2) + pow(finalSlip[2],2));
    const double slipCoef = 2.0 * _KinSrcRamp::maxAcc(finalSlipMag, riseTime, impulseDuration) / impulseDuration;

    double slipAccTimeFn = 0.0;
    if (t-t0 < 0.0) {
        slipAccTimeFn = 0.0;
    } else if (t-t0 <= 0.5*impulseDuration) {
        slipAccTimeFn = t;
    } else if (t-t0 <= impulseDuration) {
        slipAccTimeFn = impulseDuration-t;
    } else if (t-t0 <= riseTime-impulseDuration) {
        slipAccTimeFn = 0.0;
    } else if (t-t0 <= riseTime-0.5*impulseDuration) {
        slipAccTimeFn = -(t-(riseTime-impulseDuration));
    } else if (t-t0 <= riseTime) {
        slipAccTimeFn = t-riseTime;
    } else {
        slipAccTimeFn = 0.0;
    } // if/else

    for (PylithInt i = 0; i < dim; ++i) {
        slipAcc[i] = slipCoef * slipAccTimeFn * finalSlip[i] / finalSlipMag;
    } // for
} // slipAccFn


// ------------------------------------------------------------------------------------------------
// Preinitialize earthquake source. Set names/sizes of auxiliary subfields.
void
pylith::faults::KinSrcRamp::_auxiliaryFieldSetup(const spatialdata::units::Nondimensional& normalizer,
                                                 const spatialdata::geocoords::CoordSys* cs) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("_auxiliaryFieldSetup()");

    assert(_auxiliaryFactory);
    assert(cs);
    _auxiliaryFactory->initialize(_auxiliaryField, normalizer, cs->getSpaceDim());

    // :ATTENTION: The order for adding subfields must match the order of the auxiliary fields in the slip time
    // function.

    _auxiliaryFactory->addInitiationTime(); // 0
    _auxiliaryFactory->addFinalSlip(); // 1
    _auxiliaryFactory->addRiseTime(); // 2
    _auxiliaryFactory->addImpulseDuration(); // 3

    _slipFnKernel = pylith::faults::KinSrcRamp::slipFn;
    _slipRateFnKernel = pylith::faults::KinSrcRamp::slipRateFn;
    _slipAccFnKernel = pylith::faults::KinSrcRamp::slipAccFn;

    PYLITH_METHOD_END;
} // _auxiliaryFieldSetup


// ------------------------------------------------------------------------------------------------
// Maximum value of acceleration impulse in smoothed ramp slip time function.
double
pylith::faults::_KinSrcRamp::maxAcc(const double finalSlip,
                                    const double riseTime,
                                    const double impulseDuration) {
    const double a =
        1.0/6.0*pow(riseTime,3)
        -1.0/3.0*pow(riseTime-0.5*impulseDuration,3)
        +0.5*riseTime*pow(riseTime-0.5*impulseDuration,2)
        -0.5*pow(riseTime,2)*(riseTime-0.5*impulseDuration)
        +0.5*(riseTime-impulseDuration)*pow(riseTime-0.5*impulseDuration,2)
        -0.5*pow(riseTime-impulseDuration,2)*(riseTime-0.5*impulseDuration)
        +0.25*pow(impulseDuration,2)*(riseTime-0.5*impulseDuration)
        +1.0/6.0*pow(riseTime-impulseDuration,3)
        -0.125*pow(impulseDuration,3);

    return 0.5 * finalSlip * impulseDuration / a;
} // maxAcc


// End of file
