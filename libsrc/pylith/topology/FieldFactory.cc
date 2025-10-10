// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================

#include <portinfo>

#include "pylith/topology/FieldFactory.hh" // implementation of object methods

#include "pylith/topology/Field.hh" // HOLDSA AuxiliaryField

#include "pylith/scales/Scales.hh" // USES Scales

#include "pylith/utils/error.hh" // USES PYLITH_METHOD*
#include "pylith/utils/journals.hh" // USES PYLITH_JOURNAL*

#include <cassert>

// ---------------------------------------------------------------------------------------------------------------------
// Default constructor.
pylith::topology::FieldFactory::FieldFactory(void) :
    _field(NULL),
    _defaultDescription(NULL),
    _scales(new pylith::scales::Scales),
    _spaceDim(0) {
    GenericComponent::setName("auxiliaryfactory");
    _subfieldDiscretizations["default"] = pylith::topology::FieldBase::Discretization();
} // constructor


// ---------------------------------------------------------------------------------------------------------------------
// Destructor.
pylith::topology::FieldFactory::~FieldFactory(void) {
    _field = NULL; // :TODO: use shared pointer

    delete _defaultDescription;_defaultDescription = NULL;
    delete _scales;_scales = NULL;
} // destructor


// ---------------------------------------------------------------------------------------------------------------------
// Get number of subfield discretizations.
int
pylith::topology::FieldFactory::getNumSubfields(void) const {
    return _subfieldDiscretizations.size();
} // getNumSubfields


// ---------------------------------------------------------------------------------------------------------------------
// Set discretization information for auxiliary subfield.
void
pylith::topology::FieldFactory::setSubfieldDiscretization(const char* subfieldName,
                                                          const int basisOrder,
                                                          const int quadOrder,
                                                          const int dimension,
                                                          const bool isFaultOnly,
                                                          const pylith::topology::FieldBase::CellBasis cellBasis,
                                                          const pylith::topology::FieldBase::SpaceEnum feSpace,
                                                          const bool isBasisContinuous) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("setSubfieldDiscretization(subfieldName="<<subfieldName<<", basisOrder="<<basisOrder<<", quadOrder="<<quadOrder<<", dimension="<<dimension<<", cellBasis="<<cellBasis<<", isBasisContinuous="<<isBasisContinuous<<")");
    assert(dimension != 0);

    pylith::topology::FieldBase::Discretization feInfo;
    feInfo.basisOrder = basisOrder;
    feInfo.quadOrder = quadOrder;
    feInfo.dimension = dimension;
    feInfo.cellBasis = cellBasis;
    feInfo.isFaultOnly = isFaultOnly;
    feInfo.feSpace = feSpace;
    feInfo.isBasisContinuous = isBasisContinuous;
    _subfieldDiscretizations[subfieldName] = feInfo;

    PYLITH_METHOD_END;
} // setSubfieldDiscretization


// ---------------------------------------------------------------------------------------------------------------------
// Get discretization information for subfield.
const pylith::topology::FieldBase::Discretization&
pylith::topology::FieldFactory::getSubfieldDiscretization(const char* subfieldName) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("getSubfieldDiscretization(subfieldName="<<subfieldName<<")");

    pylith::topology::FieldBase::discretizations_map::const_iterator iter = _subfieldDiscretizations.find(subfieldName);
    if (iter != _subfieldDiscretizations.end()) {
        PYLITH_METHOD_RETURN(iter->second);
    } else { // not found so try default
        iter = _subfieldDiscretizations.find("default");
        if (iter == _subfieldDiscretizations.end()) {
            throw std::logic_error("Default discretization not set in field factory.");
        } // if
    } // if/else

    PYLITH_METHOD_RETURN(iter->second); // default
} // getSubfieldDiscretization


// ---------------------------------------------------------------------------------------------------------------------
// Initialie factory for setting up auxiliary subfields.
void
pylith::topology::FieldFactory::initialize(pylith::topology::Field* field,
                                           const pylith::scales::Scales& scales,
                                           const int spaceDim,
                                           const pylith::topology::FieldBase::Description* defaultDescription) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("initialize(field="<<field<<", scales="<<&scales<<", spaceDim="<<spaceDim<<", defaultDescription="<<defaultDescription<<")");

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
    assert(_scales);
    *_scales = scales;
    _spaceDim = spaceDim;

    assert(1 <= _spaceDim && _spaceDim <= 3);

    PYLITH_METHOD_END;
} // initialize


// End of file
