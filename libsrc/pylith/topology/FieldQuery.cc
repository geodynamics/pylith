// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2024, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================

#include <portinfo>

#include "pylith/topology/FieldQuery.hh" // implementation of class methods

#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/Mesh.hh" // USES Mesh

#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys
#include "spatialdata/spatialdb/SpatialDB.hh" // USES SpatialDB

#include "pylith/utils/error.hh" // USES PYLITH_CHECK_ERROR
#include "pylith/utils/EventLogger.hh" // USES EventLogger

namespace pylith {
    namespace topology {
        class _FieldQuery {
public:

            /** Find indices of spatial database values to use for subfield.
             *
             * Allocate buffer for query values.
             */
            static
            void findQueryIndices(FieldQuery::DBQueryContext* context,
                                  const pylith::string_vector& valuesForSubfield);

            // Logging events
            class Events {
public:

                static
                void init(void);

                static pylith::utils::EventLogger logger;
                static PylithInt queryDB;
                static PylithInt queryDBLabel;
                static PylithInt openDB;
            };

        }; // _FieldQuery
    } // topology
} // pylith
pylith::utils::EventLogger pylith::topology::_FieldQuery::Events::logger;
PylithInt pylith::topology::_FieldQuery::Events::queryDB;
PylithInt pylith::topology::_FieldQuery::Events::queryDBLabel;
PylithInt pylith::topology::_FieldQuery::Events::openDB;

// ------------------------------------------------------------------------------------------------
const PylithReal pylith::topology::FieldQuery::SCALE_TOLERANCE = 25.0;

// ------------------------------------------------------------------------------------------------
void
pylith::topology::_FieldQuery::Events::init(void) {
    logger.setClassName("FieldQuery");
    logger.initialize();
    openDB = logger.registerEvent("PL:FieldQuery:openDB");
    queryDB = logger.registerEvent("PL:FieldQuery:queryDB");
    queryDBLabel = logger.registerEvent("PL:FieldQuery:queryDBLabel");
} // init


// ----------------------------------------------------------------------
// Default constructor.
pylith::topology::FieldQuery::FieldQuery(const Field& field) :
    _field(field),
    _functions(NULL),
    _contexts(NULL),
    _contextPtrs(NULL) {
    _FieldQuery::Events::init();
} // constructor


// ----------------------------------------------------------------------
// Destructor.
pylith::topology::FieldQuery::~FieldQuery(void) {
    deallocate();
} // destructor


// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::topology::FieldQuery::deallocate(void) {
    PYLITH_METHOD_BEGIN;

    delete[] _functions;_functions = NULL;
    delete[] _contexts;_contexts = NULL;
    delete[] _contextPtrs;_contextPtrs = NULL;

    _subfieldQueries.clear();

    PYLITH_METHOD_END;
} // deallocate


// ----------------------------------------------------------------------
// Set query function information for subfield.
void
pylith::topology::FieldQuery::setQuery(const char* subfield,
                                       const char* queryValues[],
                                       const size_t numValues,
                                       convertfn_type converter,
                                       spatialdata::spatialdb::SpatialDB* db) {
    PYLITH_METHOD_BEGIN;

    assert(subfield);

    SubfieldQuery query;

    if (queryValues && (numValues > 0)) {
        query.queryValues.resize(numValues);
        for (size_t i = 0; i < numValues; ++i) {
            query.queryValues[i] = queryValues[i];
        } // for
    } else {
        const Field::SubfieldInfo& info = _field.getSubfieldInfo(subfield);
        query.queryValues = info.description.componentNames;
    } // if/else

    query.converter = converter;
    query.db = db;

    _subfieldQueries[subfield] = query;

    PYLITH_METHOD_END;
} // setQuery


// ----------------------------------------------------------------------
// Initialize query with default query information.
void
pylith::topology::FieldQuery::initializeWithDefaultQueries(void) {
    PYLITH_METHOD_BEGIN;

    _subfieldQueries.clear();

    const Field::subfields_type& subfields = _field._subfields;
    for (Field::subfields_type::const_iterator iter = subfields.begin(); iter != subfields.end(); ++iter) {
        const std::string& name = iter->first;
        const pylith::topology::Field::Description& description = iter->second.description;
        SubfieldQuery query;
        query.queryValues = description.componentNames;
        _subfieldQueries[name] = query;
    } // for

    PYLITH_METHOD_END;
} // initializeWithDefaultQueries


// ----------------------------------------------------------------------
// Query spatial database to set values in field.
void
pylith::topology::FieldQuery::openDB(spatialdata::spatialdb::SpatialDB* db,
                                     const PylithReal lengthScale) {
    PYLITH_METHOD_BEGIN;
    _FieldQuery::Events::logger.eventBegin(_FieldQuery::Events::openDB);

    // Open spatial database.
    if (db) {
        db->open();
    } // if

    // Create contexts. Need to put contexts into an array of
    // pointers, since Petsc function doesn't know the size of the
    // context.
    const Field::subfields_type& subfields = _field._subfields;
    const size_t size = subfields.size();
    delete[] _functions;_functions = (size > 0) ? new queryfn_type[size] : NULL;
    delete[] _contexts;_contexts = (size > 0) ? new DBQueryContext[size] : NULL;
    delete[] _contextPtrs;_contextPtrs = (size > 0) ? new DBQueryContext*[size] : NULL;

    for (Field::subfields_type::const_iterator iter = subfields.begin(); iter != subfields.end(); ++iter) {
        const std::string& name = iter->first;
        const PylithInt index = iter->second.index;
        assert(size_t(index) < subfields.size());
        const subfieldquery_map_type::const_iterator& query = _subfieldQueries.find(name);
        if (query != _subfieldQueries.end()) {
            spatialdata::spatialdb::SpatialDB* dbSubfield = (query->second.db) ? query->second.db : db;
            _functions[index] = (dbSubfield) ? queryDBPointFn : NULL;

            _contexts[index].converter = query->second.converter;
            _contexts[index].db = dbSubfield;
            _contexts[index].cs = _field.getMesh().getCoordSys();

            if (dbSubfield) {
                _FieldQuery::findQueryIndices(&_contexts[index], query->second.queryValues);
            } // if
        } else {
            _functions[index] = NULL;
        } // if/else
        _contexts[index].lengthScale = lengthScale;

        const pylith::topology::Field::Description& description = iter->second.description;
        _contexts[index].description = description.label;
        _contexts[index].valueScale = description.scale;
        _contexts[index].validator = description.validator;
        _contexts[index].validatorTolerance = description.validatorTolerance;

        _contextPtrs[index] = &_contexts[index];
    } // for

    _FieldQuery::Events::logger.eventEnd(_FieldQuery::Events::openDB);
    PYLITH_METHOD_END;
} // openDB


// ----------------------------------------------------------------------
// Query spatial database to set values in field.
void
pylith::topology::FieldQuery::queryDB(void) {
    PYLITH_METHOD_BEGIN;
    _FieldQuery::Events::logger.eventBegin(_FieldQuery::Events::queryDB);

    PetscErrorCode err = 0;
    PetscReal dummyTime = 0.0;
    err = DMProjectFunctionLocal(_field.getDM(), dummyTime, _functions, (void**)_contextPtrs, INSERT_ALL_VALUES,
                                 _field.getLocalVector());PYLITH_CHECK_ERROR(err);

    _FieldQuery::Events::logger.eventEnd(_FieldQuery::Events::queryDB);
    PYLITH_METHOD_END;
} // queryDB


// ----------------------------------------------------------------------
// Query spatial database to set values in field.
void
pylith::topology::FieldQuery::queryDBLabel(const char* labelName,
                                           const PylithInt labelValue) {
    PYLITH_METHOD_BEGIN;
    _FieldQuery::Events::logger.eventBegin(_FieldQuery::Events::queryDBLabel);

    PetscErrorCode err = 0;
    PetscReal dummyTime = 0.0;

    const Field::subfields_type& subfields = _field._subfields;
    const size_t numSubfields = subfields.size();
    pylith::int_array subfieldIndices(numSubfields);
    size_t i = 0;
    for (Field::subfields_type::const_iterator iter = subfields.begin(); iter != subfields.end(); ++iter, ++i) {
        subfieldIndices[i] = iter->second.index;
    } // for

    PetscDMLabel dmLabel = NULL;
    err = DMGetLabel(_field.getDM(), labelName, &dmLabel);PYLITH_CHECK_ERROR(err);assert(dmLabel);
    err = DMProjectFunctionLabelLocal(_field.getDM(), dummyTime, dmLabel, 1, &labelValue,
                                      numSubfields, &subfieldIndices[0], _functions, (void**)_contextPtrs,
                                      INSERT_ALL_VALUES, _field.getLocalVector());PYLITH_CHECK_ERROR(err);

    _FieldQuery::Events::logger.eventEnd(_FieldQuery::Events::queryDBLabel);
    PYLITH_METHOD_END;
} // queryDBLabel


// ----------------------------------------------------------------------
// Query spatial database to set values in field.
void
pylith::topology::FieldQuery::closeDB(spatialdata::spatialdb::SpatialDB* db) {
    PYLITH_METHOD_BEGIN;

    delete[] _functions;_functions = NULL;
    delete[] _contexts;_contexts = NULL;
    delete[] _contextPtrs;_contextPtrs = NULL;

    // Close spatial database.
    if (db) {
        db->close();
    } // if

    PYLITH_METHOD_END;
} // queryDB


// ----------------------------------------------------------------------
// Generic query of values from spatial database.
PetscErrorCode
pylith::topology::FieldQuery::queryDBPointFn(PylithInt dim,
                                             PylithReal t,
                                             const PylithReal x[],
                                             PylithInt nvalues,
                                             PylithScalar* values,
                                             void* context) {
    PYLITH_METHOD_BEGIN;

    assert(x);
    assert(values);
    assert(context);

    pylith::topology::FieldQuery::DBQueryContext* queryctx = (pylith::topology::FieldQuery::DBQueryContext*)context;assert(queryctx);
    if (!queryctx->db) {
        PYLITH_METHOD_RETURN(0);
    } // if

    // Dimensionalize query location coordinates.
    assert(queryctx->lengthScale > 0);
    double xDim[3];
    for (int i = 0; i < dim; ++i) {
        xDim[i] = x[i] * queryctx->lengthScale;
    } // for

    // Query database.
    assert(queryctx->cs);
    const int err = queryctx->db->query(&queryctx->queryValues[0], queryctx->queryValues.size(), xDim, dim, queryctx->cs);
    if (err) {
        std::ostringstream msg;
        msg << "Could not find values for " << queryctx->description << " at (";
        for (int i = 0; i < dim; ++i) {
            msg << "  " << xDim[i];
        }
        msg << ") in spatial database '" << queryctx->db->getDescription() << "'.";
        PYLITH_ERROR_RETURN(PETSC_COMM_SELF, PETSC_ERR_LIB, msg.str().c_str());
    } // if

    // Convert database values to subfield values if converter function specified.
    if (queryctx->converter) {
        const std::string& invalidMsg = queryctx->converter(values, nvalues, queryctx->queryValues, queryctx->queryIndices);
        if (invalidMsg.length() > 0) {
            std::ostringstream msg;
            msg << "Error converting spatial database values for " << queryctx->description << " at (";
            for (int i = 0; i < dim; ++i) {
                msg << "  " << xDim[i];
            }
            msg << ") in spatial database '" << queryctx->db->getDescription() << "'. "
                << invalidMsg;
            PYLITH_ERROR_RETURN(PETSC_COMM_SELF, PETSC_ERR_LIB, msg.str().c_str());
        }
    } else {
        for (PylithInt i = 0; i < nvalues; ++i) {
            values[i] = queryctx->queryValues[queryctx->queryIndices[i]];
        } // for
    } // if/else

    // Validate subfield values if validator function was specified.
    if (queryctx->validator) {
        for (PylithInt i = 0; i < nvalues; ++i) {
            const std::string& invalidMsg = queryctx->validator(values[i], queryctx->valueScale, queryctx->validatorTolerance);
            if (invalidMsg.length() > 0) {
                std::ostringstream msg;
                msg << "Found invalid value (" << values[i] << ") for " << queryctx->description
                    << " at location (";
                for (int i = 0; i < dim; ++i) {
                    msg << "  " << xDim[i];
                } // for
                msg << ") from spatial database '" << queryctx->db->getDescription() << "'. ";
                msg << invalidMsg;
                PYLITH_ERROR_RETURN(PETSC_COMM_SELF, PETSC_ERR_LIB, msg.str().c_str());
            } // if
        } // for
    } // if

    // Nondimensionalize values
    assert(queryctx->valueScale > 0);
    for (int i = 0; i < nvalues; ++i) {
        values[i] /= queryctx->valueScale;
    } // for

    PYLITH_METHOD_RETURN(0);
} // queryDBPointFn


// ----------------------------------------------------------------------
std::string
pylith::topology::FieldQuery::validatorPositive(const PylithReal value,
                                                const PylithReal scale,
                                                const PylithReal tolerance) {
    std::string errorMsg;
    if (value <= 0.0) {
        errorMsg = std::string("Value must be positive.");
    } else if ((scale > 0.0) && (tolerance > 0.0)) {
        const PylithReal minValue = scale / tolerance;
        const PylithReal maxValue = scale * tolerance;
        if ((value < minValue) || (value > maxValue)) {
            std::ostringstream msg;
            msg << "Value outside range [" << minValue << ", " << maxValue << "] for nondimensionalization. "
                << "You likely need to adjust the scales for nondimensionalization.";
            errorMsg = msg.str();
        } // if
    } // if
    return errorMsg;
} // validatorPositive


// ----------------------------------------------------------------------
std::string
pylith::topology::FieldQuery::validatorNonnegative(const PylithReal value,
                                                   const PylithReal scale,
                                                   const PylithReal tolerance) {
    std::string errorMsg;
    if (value < 0.0) {
        errorMsg = std::string("Value must be non-negative.");
    } else if ((value > 0) && (scale > 0.0) && (tolerance > 0)) {
        const PylithReal minValue = scale / tolerance;
        const PylithReal maxValue = scale * tolerance;
        if ((value < minValue) || (value > maxValue)) {
            std::ostringstream msg;
            msg << "Value outside range [" << minValue << ", " << maxValue << "] for nondimensionalization. "
                << "You likely need to adjust the scales for nondimensionalization.";
            errorMsg = msg.str();
        } // if
    } // if/else

    return errorMsg;
} // validatorNonnegative


// ----------------------------------------------------------------------
std::string
pylith::topology::FieldQuery::validatorScale(const PylithReal value,
                                             const PylithReal scale,
                                             const PylithReal tolerance) {
    std::string errorMsg;
    if ((scale > 0.0) && (tolerance > 0)) {
        const PylithReal minValue = scale / tolerance;
        const PylithReal maxValue = scale * tolerance;
        if ((fabs(value) < minValue) || (fabs(value) > maxValue)) {
            std::ostringstream msg;
            msg << "Absolute value outside range [" << minValue << ", " << maxValue << "] for nondimensionalization. "
                << "You likely need to adjust the scales for nondimensionalization.";
            errorMsg = msg.str();
        } // if
    } // if/else

    return errorMsg;
} // validatorNonnegative


// ----------------------------------------------------------------------
// Find indices of spatial database values to use for subfield. Allocate buffer for query values.
void
pylith::topology::_FieldQuery::findQueryIndices(FieldQuery::DBQueryContext* context,
                                                const pylith::string_vector& valuesForSubfield) {
    assert(context);

    const char** dbValues = NULL;
    size_t numDBValues = 0;
    context->db->getNamesDBValues(&dbValues, &numDBValues);

    const size_t numValues = valuesForSubfield.size();
    context->queryIndices.resize(numValues);
    for (size_t iValue = 0; iValue < numValues; ++iValue) {
        bool foundName = false;
        for (size_t index = 0; index < numDBValues; ++index) {
            if (0 == strcasecmp(dbValues[index], valuesForSubfield[iValue].c_str())) {
                foundName = true;
                context->queryIndices[iValue] = index;
                break;
            } // if
        } // for
        if (!foundName) {
            std::ostringstream msg;
            if (0 == numDBValues) {
                delete dbValues;dbValues = NULL;
                msg << "No values found in spatial database '"
                    << context->db->getDescription() << "'. Did you forget to open the database?";
                throw std::logic_error(msg.str());
            } // if
            msg << "Could not find value '" << valuesForSubfield[iValue] << "' in spatial database '"
                << context->db->getDescription() << "'. Available values are:";
            for (size_t iValueDB = 0; iValueDB < numDBValues; ++iValueDB) {
                msg << "\n  " << dbValues[iValueDB];
            } // for
            msg << "\n";
            delete dbValues;dbValues = NULL;
            throw std::out_of_range(msg.str());
        } // if
    } // for
    delete[] dbValues;dbValues = NULL;

    context->queryValues.resize(numDBValues);
} // findQueryIndices


// End of file
