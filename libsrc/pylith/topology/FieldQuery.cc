// -*- C++ -*-
//
// ======================================================================
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
// ======================================================================
//

#include <portinfo>

#include "FieldQuery.hh" // implementation of class methods

#include "Field.hh" // USES Field
#include "Mesh.hh" // USES Mesh

#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys
#include "spatialdata/spatialdb/SpatialDB.hh" // USES SpatialDB

#include "pylith/utils/error.hh" // USES PYLITH_CHECK_ERROR

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

        }; // _FieldQuery
    } // topology
} // pylith

// ----------------------------------------------------------------------
// Default constructor.
pylith::topology::FieldQuery::FieldQuery(const Field& field) :
    _field(field),
    _functions(NULL),
    _contexts(NULL),
    _contextPtrs(NULL) {}


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
                                       convertfn_type convertFn,
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
        const Field::SubfieldInfo& info = _field.subfieldInfo(subfield);
        query.queryValues = info.description.componentNames;
    } // if/else

    query.convertFn = convertFn;
    query.db = db;

    _subfieldQueries[subfield] = query;

    PYLITH_METHOD_END;
} // setQuery


#if 0
// ----------------------------------------------------------------------
// Get query function information for subfield.
const pylith::topology::FieldQuery::queryfn_type
pylith::topology::FieldQuery::queryFn(const char* subfield) const {
    PYLITH_METHOD_BEGIN;

    assert(subfield);
    const queryfn_map_type::const_iterator& iter = _subfieldQueries.find(subfield);
    if (iter == _subfieldQueries.end()) {
        std::ostringstream msg;
        msg << "Could not find query function for subfield '" << subfield << "'." << std::endl;
        throw std::logic_error(msg.str());
    } // if

    PYLITH_METHOD_RETURN(*iter->second);
} // queryFn


#endif

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


#if 0
// ----------------------------------------------------------------------
// Get spatial database used to get values for subfield.
const spatialdata::spatialdb::SpatialDB*
pylith::topology::FieldQuery::queryDB(const char* subfield) const {
    PYLITH_METHOD_BEGIN;

    assert(subfield);
    const querydb_map_type::const_iterator& iter = _queryDBs.find(subfield);
    if (iter == _queryDBs.end()) {
        std::ostringstream msg;
        msg << "Could not find spatial database for subfield '" << subfield << "'." << std::endl;
        throw std::logic_error(msg.str());
    } // if

    PYLITH_METHOD_RETURN(iter->second);
} // queryDB


#endif

#if 0
// ----------------------------------------------------------------------
// Get array of query functions.
pylith::topology::FieldQuery::queryfn_type*
pylith::topology::FieldQuery::functions(void) const {
    return _functions;
} // functions


// ----------------------------------------------------------------------
// Get array of pointers to contexts.
const pylith::topology::FieldQuery::DBQueryContext* const*
pylith::topology::FieldQuery::contextPtrs(void) const {
    return _contextPtrs;
} // contextPtrs


#endif

// ----------------------------------------------------------------------
// Query spatial database to set values in field.
void
pylith::topology::FieldQuery::openDB(spatialdata::spatialdb::SpatialDB* db,
                                     const PylithReal lengthScale) {
    PYLITH_METHOD_BEGIN;

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

            _contexts[index].db = dbSubfield;
            _contexts[index].cs = _field.mesh().getCoordSys();

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

        _contextPtrs[index] = &_contexts[index];
    } // for

    // Open spatial database.
    if (db) {
        db->open();
    } // if

    PYLITH_METHOD_END;
} // openDB


// ----------------------------------------------------------------------
// Query spatial database to set values in field.
void
pylith::topology::FieldQuery::queryDB(void) {
    PYLITH_METHOD_BEGIN;

    PetscErrorCode err = 0;
    PetscReal dummyTime = 0.0;
    err = DMProjectFunctionLocal(_field.dmMesh(), dummyTime, _functions, (void**)_contextPtrs, INSERT_ALL_VALUES,
                                 _field.localVector());PYLITH_CHECK_ERROR(err);

    PYLITH_METHOD_END;
} // queryDB


// ----------------------------------------------------------------------
// Query spatial database to set values in field.
void
pylith::topology::FieldQuery::queryDBLabel(const char* labelName,
                                           const PylithInt labelValue) {
    PYLITH_METHOD_BEGIN;

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
    err = DMGetLabel(_field.dmMesh(), labelName, &dmLabel);PYLITH_CHECK_ERROR(err);
    err = DMProjectFunctionLabelLocal(_field.dmMesh(), dummyTime, dmLabel, 1, &labelValue,
                                      numSubfields, &subfieldIndices[0], _functions, (void**)_contextPtrs,
                                      INSERT_ALL_VALUES, _field.localVector());PYLITH_CHECK_ERROR(err);

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
        msg << ") in spatial database '" << queryctx->db->getLabel() << "'.";
        PYLITH_SET_ERROR(PETSC_COMM_SELF, PETSC_ERR_LIB, msg.str().c_str());
        PYLITH_METHOD_RETURN(1);
    } // if

    // Convert database values to subfield values if converter function specified.
    assert(queryctx->queryIndices.size() == nvalues);
    if (queryctx->converter) {
        PetscErrorCode err = queryctx->converter(values, nvalues, queryctx->queryValues, queryctx->queryIndices);
        if (err) { PYLITH_METHOD_RETURN(1); }
    } else {
        for (PylithInt i = 0; i < nvalues; ++i) {
            values[i] = queryctx->queryValues[queryctx->queryIndices[i]];
        } // for
    } // if/else

    // Validate subfield values if validator function was specified.
    if (queryctx->validator) {
        for (PylithInt i = 0; i < nvalues; ++i) {
            const char* invalidMsg = queryctx->validator(values[i]);
            if (invalidMsg) {
                std::ostringstream msg;
                msg << "Found invalid value for " << queryctx->description << " (" << values[i] << ") at location (";
                for (int i = 0; i < dim; ++i) {
                    msg << "  " << xDim[i];
                }
                msg << ") from spatial database '" << queryctx->db->getLabel() << "'. ";
                msg << invalidMsg;
                PYLITH_SET_ERROR(PETSC_COMM_SELF, PETSC_ERR_LIB, msg.str().c_str());
                PYLITH_METHOD_RETURN(1);
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
const char*
pylith::topology::FieldQuery::validatorPositive(const PylithReal value) {
    return (value > 0.0) ? NULL : "Value must be positive.";
} // validatorPositive


// ----------------------------------------------------------------------
const char*
pylith::topology::FieldQuery::validatorNonnegative(const PylithReal value) {
    return (value >= 0.0) ? NULL : "Value must be nonnegative.";
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
            msg << "Could not find value '" << valuesForSubfield[iValue] << "' in spatial database '"
                << context->db->getLabel() << "'. Available values are:";
            for (size_t iValueDB = 0; iValueDB < numDBValues; ++iValueDB) {
                msg << "\n  " << dbValues[iValueDB];
            } // for
            msg << "\n";
            throw std::out_of_range(msg.str());
        } // if
    } // for

    context->queryValues.resize(numDBValues);
} // findQueryIndices


// End of file
