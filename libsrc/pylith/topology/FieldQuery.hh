// -*- C++ -*-
//
// ======================================================================
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
// ======================================================================
//

/**
 * @file libsrc/topology/FieldQuery.hh
 *
 * @brief Set field via query of spatial database, etc.
 */

#if !defined(pylith_topology_fieldquery_hh)
#define pylith_topology_fieldquery_hh

#include "pylith/topology/topologyfwd.hh" // forward declarations

#include "pylith/topology/FieldBase.hh" // HASA validatorfn_type
#include "pylith/testing/testingfwd.hh" // USES FieldTester
#include "pylith/utils/utilsfwd.hh" // USES EventLogger

#include "spatialdata/spatialdb/spatialdbfwd.hh" // HOLDSA SpatialDB
#include "spatialdata/geocoords/geocoordsfwd.hh" // USES CoordSys

#include <map> // HOLDSA std::map
#include <string> // HASA std::string

namespace pylith {
    namespace feassemble {
        class TestAuxiliaryFactory;
    } // feassemble
} // pylith

class pylith::topology::FieldQuery {
    friend class _FieldQuery;
    friend class TestFieldQuery; // unit testing
    friend class pylith::testing::FieldTester; // unit testing
    friend class pylith::feassemble::TestAuxiliaryFactory; // unit testing

    // PUBLIC TYPEDEF ///////////////////////////////////////////////////////
public:

    /** Function prototype for converter functions.
     *
     * @param[out] values Values for subfield.
     * @param[in] nvalues Number of values for subfield.
     * @param[in] Array of values from spatial database query.
     * @param[in] Indices of values from spatial database to use for computing subfield values.
     */
    typedef std::string (*convertfn_type)(PylithScalar[],
                                          const PylithInt,
                                          const pylith::scalar_array,
                                          const pylith::int_array);

    // PUBLIC MEMBERS ///////////////////////////////////////////////////////
public:

    /** Default constructor.
     *
     * @param field Field associated with query.
     */
    FieldQuery(const Field& field);

    /// Destructor.
    ~FieldQuery(void);

    /// Deallocate PETSc and local data structures.
    void deallocate(void);

    /** Set query information for subfield.
     *
     * The default is to use the database passed in the call to openDB().
     * Passing in a spatial database in this function overrides use of the
     * default database. If the names of the values to query in the spatial database are not given via queryValues,
     * then the names of the subfield components are used.
     *
     * @param[in] subfield Name of subfield.
     * @param[in] queryValues Array of names of spatial database values for subfield.
     * @param[in] numValues Size of names array.
     * @param[in] converter Function to convert spatial database values to subfield value (optional).
     * @param[in] db Spatial database to query (optional).
     */
    void setQuery(const char* subfield,
                  const char* queryValues[]=NULL,
                  const size_t numValues=0,
                  convertfn_type converter=NULL,
                  spatialdata::spatialdb::SpatialDB* db=NULL);

    /// Initialize query with default query information.
    void initializeWithDefaultQueries(void);

    /** Open spatial database query for setting values in field.
     *
     * @param db Spatial database to query.
     * @param lengthScale Length scale for dimensionalization of
     * location coordinates.
     */
    void openDB(spatialdata::spatialdb::SpatialDB* db,
                const PylithReal lengthScale);

    /// Query spatial database to set values in field.
    void queryDB(void);

    /** Query spatial database for points in label to set values in field.
     *
     * @param[in] name Name of label.
     * @param[in] value Value of label.
     */
    void queryDBLabel(const char* name,
                      const PylithInt value=1);

    /** Close spatial database query for setting values in field.
     *
     * @param db Spatial database to query.
     */
    void closeDB(spatialdata::spatialdb::SpatialDB* db);

    /** Query of values from spatial databaseat point.
     *
     * Includes nondimensionalization but no conversion of values.
     *
     * @param[in] dim Spatial dimension.
     * @param[in] t Current time.
     * @param[in] x Coordinates (nondimensioned) of point location for query.
     * @param[in] nvalues Size of values array.
     * @param[out] values Array of values to be returned.
     * @param[in] context Query context.
     * @returns PETSc error code (0 for success).
     */
    static
    PetscErrorCode queryDBPointFn(PylithInt dim,
                                  PylithReal t,
                                  const PylithReal x[],
                                  PylithInt nvalues,
                                  PylithScalar* values,
                                  void* context);

    /** Validator for positive values.
     *
     * @param[in] value Value to validate.
     * @returns Error message if not positive, NULL otherwise.
     */
    static
    const char* validatorPositive(const PylithReal value);

    /** Validator for nonnegative values.
     *
     * @param[in] value Value to validate.
     * @returns Error message if negative, NULL otherwise.
     */
    static
    const char* validatorNonnegative(const PylithReal value);

    // PRIVATE TYPEDEFS /////////////////////////////////////////////////////
private:

    struct SubfieldQuery {
        pylith::string_vector queryValues; ///< Values to use from spatial database.
        convertfn_type converter; ///< Function to convert spatial database values to subfield values.
        spatialdata::spatialdb::SpatialDB* db; ///< Spatial database to query.

        SubfieldQuery(void) :
            converter(NULL),
            db(NULL) {}


    }; // SubfieldQuery

    typedef std::map<std::string, SubfieldQuery> subfieldquery_map_type;

    /// Function prototype for query functions.
    typedef PetscErrorCode (*queryfn_type)(PylithInt,
                                           PylithReal,
                                           const PylithReal[],
                                           PylithInt,
                                           PylithScalar*,
                                           void*);

    /// Context for spatial database queries.
    struct DBQueryContext {
        spatialdata::spatialdb::SpatialDB* db; ///< Spatial database.
        const spatialdata::geocoords::CoordSys* cs; ///< Coordinate system of input point locations.
        PylithReal lengthScale; ///< Length scale for dimensionalizing coordinates.
        PylithReal valueScale; ///< Scale for dimensionalizing values for subfield.
        std::string description; ///< Name of value;
        pylith::scalar_array queryValues; ///< Values returned by spatial database query;
        pylith::int_array queryIndices; ///< Indices of spatial database values to use for subfield.
        convertfn_type converter; ///< Function to convert values to subfield (optional).
        pylith::topology::FieldBase::validatorfn_type validator; ///< Function to validate values (optional).
        pylith::utils::EventLogger* logger;

        DBQueryContext(void) :
            db(NULL),
            cs(NULL),
            lengthScale(1.0),
            valueScale(1.0),
            description("unknown"),
            converter(NULL),
            validator(NULL),
            logger(NULL) {}


    }; // DBQueryStruct

    // PRIVATE MEMBERS //////////////////////////////////////////////////////
private:

    const pylith::topology::Field& _field; ///< Field associated with query.

    subfieldquery_map_type _subfieldQueries;

    queryfn_type* _functions; ///< Functions implementing queries.
    DBQueryContext* _contexts; ///< Contexts for performing query for each subfield.
    DBQueryContext** _contextPtrs; ///< Array of pointers to contexts.

    pylith::utils::EventLogger* _logger;

    // NOT IMPLEMENTED //////////////////////////////////////////////////////
private:

    FieldQuery(const FieldQuery&); ///< Not implemented
    const FieldQuery& operator=(const FieldQuery&); ///< Not implemented

}; // FieldQuery

#endif // pylith_topology_fieldquery_hh

// End of file
