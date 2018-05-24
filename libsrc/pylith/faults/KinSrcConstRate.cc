// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2017 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "KinSrcConstRate.hh" // implementation of object methods

#include "pylith/faults/KinSrcAuxiliaryFactory.hh" // USES KinSrcAuxiliaryFactory

#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_*
#include "pylith/utils/error.hh" // USES PYLITH_METHOD_BEGIN/END

#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys

#include <cassert> // USES assert()

// ----------------------------------------------------------------------
const char* pylith::faults::KinSrcConstRate::_pyreComponent = "kinsrcconstrate";


// ----------------------------------------------------------------------
// Default constructor.
pylith::faults::KinSrcConstRate::KinSrcConstRate(void)
{ // constructor
    pylith::utils::PyreComponent::name(_pyreComponent);
} // constructor


// ----------------------------------------------------------------------
// Destructor.
pylith::faults::KinSrcConstRate::~KinSrcConstRate(void) {}

// ----------------------------------------------------------------------
// Slip time function kernel.
void
pylith::faults::KinSrcConstRate::slipFn(const PylithInt dim,
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
    const PylithInt _numA = 2;

    assert(_numA == numA);
    assert(aOff);
    assert(a);
    assert(slip);

    const PylithInt i_slipTime = 0;
    const PylithInt i_slipRate = 1;
    const PylithScalar slipTime = a[aOff[i_slipTime]];
    const PylithScalar* slipRate = &a[aOff[i_slipRate]];

    const PylithInt i_originTime = 0;
    const PylithScalar originTime = constants[i_originTime];
    const PylithScalar t0 = originTime + slipTime;

    if (t >= t0) {
        for (PylithInt i = 0; i < dim; ++i) {
            slip[i] += slipRate[i] * (t - t0);
        } // for
    } // if

} // slipFn

// ----------------------------------------------------------------------
// Preinitialize earthquake source. Set names/sizes of auxiliary subfields.
void
pylith::faults::KinSrcConstRate::_auxFieldSetup(const spatialdata::units::Nondimensional& normalizer,
                                                const spatialdata::geocoords::CoordSys* cs) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("_auxFieldSetup()");

    assert(_auxFactory);
    assert(cs);
    _auxFactory->initialize(_auxField, normalizer, cs->spaceDim());

    // :ATTENTION: The order for adding subfields must match the order of the auxiliary fields in the slip time function kernel.

    _auxFactory->slipTime(); // 0
    _auxFactory->slipRate(); // 1

    _slipFnKernel = pylith::faults::KinSrcConstRate::slipFn;

    PYLITH_METHOD_END;
} // _auxFieldSetup


// End of file
