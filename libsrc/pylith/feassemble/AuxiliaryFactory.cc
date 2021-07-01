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

#include "AuxiliaryFactory.hh" // implementation of object methods

#include "pylith/topology/Field.hh" // HOLDSA AuxiliaryField
#include "pylith/topology/FieldQuery.hh" // USES FieldQuery

#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include "pylith/utils/error.hh" // USES PYLITH_METHOD*
#include "pylith/utils/journals.hh" // USES PYLITH_JOURNAL*

#include <cassert>

// ---------------------------------------------------------------------------------------------------------------------
// Default constructor.
pylith::feassemble::AuxiliaryFactory::AuxiliaryFactory(void) :
    _queryDB(NULL),
    _fieldQuery(NULL) {
    GenericComponent::setName("auxiliaryfactory");
} // constructor


// ---------------------------------------------------------------------------------------------------------------------
// Destructor.
pylith::feassemble::AuxiliaryFactory::~AuxiliaryFactory(void) {
    _queryDB = NULL; // :TODO: use shared pointer

    delete _fieldQuery;_fieldQuery = NULL;
} // destructor


// ---------------------------------------------------------------------------------------------------------------------
// Set database for filling auxiliary subfields.
void
pylith::feassemble::AuxiliaryFactory::setQueryDB(spatialdata::spatialdb::SpatialDB* value) {
    _queryDB = value;
} // setQueryDB


// ---------------------------------------------------------------------------------------------------------------------
// Get database for filling auxiliary subfields.
const spatialdata::spatialdb::SpatialDB*
pylith::feassemble::AuxiliaryFactory::getQueryDB(void) const {
    return _queryDB;
} // getQueryDB


// ---------------------------------------------------------------------------------------------------------------------
// Initialie factory for setting up auxiliary subfields.
void
pylith::feassemble::AuxiliaryFactory::initialize(pylith::topology::Field* field,
                                                 const spatialdata::units::Nondimensional& normalizer,
                                                 const int spaceDim,
                                                 const pylith::topology::FieldBase::Description* defaultDescription) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("initialize(field="<<field<<", normalizer="<<&normalizer<<", spaceDim="<<spaceDim<<", defaultDescription="<<defaultDescription<<")");

    FieldFactory::initialize(field, normalizer, spaceDim, defaultDescription);
    delete _fieldQuery;_fieldQuery = new pylith::topology::FieldQuery(*field);

    PYLITH_METHOD_END;
} // initialize


// ---------------------------------------------------------------------------------------------------------------------
// Initialize subfields.
void
pylith::feassemble::AuxiliaryFactory::setValuesFromDB(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("setValuesFromDB()");

    assert(_normalizer);

    if (_queryDB) {
        assert(_fieldQuery);
        _fieldQuery->openDB(_queryDB, _normalizer->getLengthScale());
        _fieldQuery->queryDB();
        _fieldQuery->closeDB(_queryDB);
    } else {
        PYLITH_JOURNAL_ERROR("Unknown case for filling auxiliary subfields.");
        throw std::logic_error("Unknown case for filling auxiliary subfields.");
    } // if/else

    delete _fieldQuery;_fieldQuery = NULL;
    _field = NULL;

    // this->view("AUXILIARY FIELDS"); // :DEBUGGING: TEMPORARY

    PYLITH_METHOD_END;
} // setValuesFromDB


// ---------------------------------------------------------------------------------------------------------------------
// Set query function for subfield.
void
pylith::feassemble::AuxiliaryFactory::setSubfieldQuery(const char* subfieldName,
                                                       const char* namesDBValues[],
                                                       const size_t numDBValues,
                                                       pylith::topology::FieldQuery::convertfn_type convertFn,
                                                       spatialdata::spatialdb::SpatialDB* db) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("setSubfieldQuery(subfieldName="<<subfieldName<<", namesDBValues="<<namesDBValues<<", numDBValues="<<numDBValues<<", convertFn="<<convertFn<<", db="<<db<<")");

    assert(_fieldQuery);
    _fieldQuery->setQuery(subfieldName, namesDBValues, numDBValues, convertFn, db);

    PYLITH_METHOD_END;
} // _setSubfieldQueryFn


// End of file
