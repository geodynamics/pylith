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

#include "AuxiliaryFactory.hh" // implementation of object methods

#include "pylith/topology/Field.hh" // HOLDSA AuxiliaryField
#include "pylith/topology/FieldQuery.hh" // USES FieldQuery

#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include "pylith/utils/journals.hh" // USES PYLITH_JOURNAL*

#include <cassert>

// ---------------------------------------------------------------------------------------------------------------------
// Default constructor.
pylith::feassemble::AuxiliaryFactory::AuxiliaryFactory(void) :
    _field(NULL),
    _defaultDescription(NULL),
    _normalizer(new spatialdata::units::Nondimensional),
    _spaceDim(0),
    _queryDB(NULL),
    _fieldQuery(NULL) { // constructor
    _subfieldDiscretizations["default"] = pylith::topology::FieldBase::Discretization();
} // constructor


// ---------------------------------------------------------------------------------------------------------------------
// Destructor.
pylith::feassemble::AuxiliaryFactory::~AuxiliaryFactory(void) {
    _field = NULL; // :TODO: use shared pointer
    _queryDB = NULL; // :TODO: use shared pointer

    delete _defaultDescription;_defaultDescription = NULL;
    delete _normalizer;_normalizer = NULL;
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
pylith::feassemble::AuxiliaryFactory::getQueryDB(void) {
    return _queryDB;
} // getQueryDB


// ---------------------------------------------------------------------------------------------------------------------
// Set discretization information for auxiliary subfield.
void
pylith::feassemble::AuxiliaryFactory::setSubfieldDiscretization(const char* subfieldName,
                                                                const int basisOrder,
                                                                const int quadOrder,
                                                                const int dimension,
                                                                const bool isBasisContinuous,
                                                                const pylith::topology::FieldBase::SpaceEnum feSpace) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("setSubfieldDiscretization(subfieldName="<<subfieldName<<", basisOrder="<<basisOrder<<", quadOrder="<<quadOrder<<", isBasisContinuous="<<isBasisContinuous<<")");

    pylith::topology::FieldBase::Discretization feInfo;
    feInfo.basisOrder = basisOrder;
    feInfo.quadOrder = quadOrder;
    feInfo.isBasisContinuous = isBasisContinuous;
    feInfo.dimension = dimension;
    feInfo.feSpace = feSpace;
    _subfieldDiscretizations[subfieldName] = feInfo;

    PYLITH_METHOD_END;
} // setSubfieldDiscretization


// ---------------------------------------------------------------------------------------------------------------------
// Get discretization information for subfield.
const pylith::topology::FieldBase::Discretization&
pylith::feassemble::AuxiliaryFactory::getSubfieldDiscretization(const char* subfieldName) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("getSubfieldDiscretization(subfieldName="<<subfieldName<<")");

    pylith::topology::FieldBase::discretizations_map::const_iterator iter = _subfieldDiscretizations.find(subfieldName);
    if (iter != _subfieldDiscretizations.end()) {
        PYLITH_METHOD_RETURN(iter->second);
    } else { // not found so try default
        iter = _subfieldDiscretizations.find("default");
        if (iter == _subfieldDiscretizations.end()) {
            throw std::logic_error("Default discretization not set for auxiliary fields.");
        } // if
    } // if/else

    PYLITH_METHOD_RETURN(iter->second); // default
} // getSubfieldDiscretization


// ---------------------------------------------------------------------------------------------------------------------
// Initialie factory for setting up auxiliary subfields.
void
pylith::feassemble::AuxiliaryFactory::initialize(pylith::topology::Field* field,
                                                 const spatialdata::units::Nondimensional& normalizer,
                                                 const int spaceDim,
                                                 const pylith::topology::FieldBase::Description* defaultDescription) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("initialize(field="<<field<<", normalizer="<<&normalizer<<", spaceDim="<<spaceDim<<", defaultDescription="<<defaultDescription<<")");

    assert(field);

    _field = field;
    if (defaultDescription) {
        if (!_defaultDescription) {
            _defaultDescription = new pylith::topology::FieldBase::Description;assert(_defaultDescription);
        } // if
        *_defaultDescription = *defaultDescription;
    } else {
        delete _defaultDescription;_defaultDescription = NULL;
    } // if/else
    assert(_normalizer);
    *_normalizer = normalizer;
    _spaceDim = spaceDim;

    assert(1 <= _spaceDim && _spaceDim <= 3);

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
        _fieldQuery->openDB(_queryDB, _normalizer->lengthScale());
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
pylith::feassemble::AuxiliaryFactory::_setSubfieldQueryFn(const char* subfieldName,
                                                          pylith::topology::FieldQuery::queryfn_type fn,
                                                          spatialdata::spatialdb::SpatialDB* db) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("_setSubfieldQueryFn(subfieldName="<<subfieldName<<", queryfn_type="<<fn<<")");

    assert(_fieldQuery);
    _fieldQuery->queryFn(subfieldName, fn, db);

    PYLITH_METHOD_END;
} // _setSubfieldQueryFn


// End of file
