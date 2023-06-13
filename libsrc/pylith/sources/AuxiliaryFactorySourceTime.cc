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
// Copyright (c) 2010-2015 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "AuxiliaryFactorySourceTime.hh" // implementation of object methods

#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/FieldQuery.hh" // HOLDSA FieldQuery

#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include "pylith/utils/error.hh" // USES PYLITH_METHOD*
#include "pylith/utils/journals.hh" // USES PYLITH_JOURNAL*

#include <cassert>

// ---------------------------------------------------------------------------------------------------------------------
// Default constructor.
pylith::sources::AuxiliaryFactorySourceTime::AuxiliaryFactorySourceTime(void) {
    GenericComponent::setName("auxiliaryfactorysourcetime");
} // constructor

// ---------------------------------------------------------------------------------------------------------------------
// Destructor.
pylith::sources::AuxiliaryFactorySourceTime::~AuxiliaryFactorySourceTime(void) {}


// ---------------------------------------------------------------------------------------------------------------------
// Add center frequency of source time function to auxiliary fields.
void
pylith::sources::AuxiliaryFactorySourceTime::addCenterFrequency(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addCenterFrequency(void)");

    const char* subfieldName = "center_frequency";

    pylith::topology::Field::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = subfieldName;
    description.scale = 1.0;
    description.validator = pylith::topology::FieldQuery::validatorNonnegative;

    _field->subfieldAdd(description, getSubfieldDiscretization(subfieldName));
    this->setSubfieldQuery(subfieldName);
    
    PYLITH_METHOD_END;
} // addRickerCenterFrequency

// End of file
