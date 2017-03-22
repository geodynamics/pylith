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

/**
 * @file libsrc/topology/FieldQuery.hh
 *
 * @brief Set field via query of spatial database, etc.
 */

#if !defined(pylith_topology_fieldquery_hh)
#define pylith_topology_fieldquery_hh

// Include directives ---------------------------------------------------
#include "pylith/topology/topologyfwd.hh" // forward declarations

#include "pylith/topology/FieldBase.hh" // HASA validatorfn_type

#include <map> // HOLDSA std::map
#include <string> // USES std::string

#include "spatialdata/spatialdb/SpatialDB.hh" // HOLDSA SpatialDB

// FieldQuery ----------------------------------------------------------------
/** @brief Set field via query of spatial database, etc.
 *
 * Each subfield is queried separately via user-defined functions to
 * allow for conversions, dimensionalization, etc.
 */
class pylith::topology::FieldQuery
{ // FieldQuery
    friend class TestFieldQuery; // unit testing

    // PUBLIC TYPEDEF ///////////////////////////////////////////////////////
public:

    /// Function prototype for query functions.
    typedef PetscErrorCode (*queryfn_type)(PylithInt, PylithReal, const PylithReal[], PylithInt, PylithScalar*, void*);

    // PUBLIC STRUCT ////////////////////////////////////////////////////////
public:

    /// Context for spatial database queries.
    struct DBQueryContext {
        spatialdata::spatialdb::SpatialDB* db; ///< Spatial database.
        const spatialdata::geocoords::CoordSys* cs; ///< Coordinate system of point locations.
        PylithReal lengthScale; ///< Length scale for dimensionalizing coordinates.
        PylithReal valueScale; ///< Scale for dimensionalizing values for subfield.
        std::string description; ///< Name of value;
        pylith::string_vector componentNames; ///< Names of components to query (optional).
        pylith::topology::FieldBase::validatorfn_type validator; ///< Function to validate values (optional).
    }; // DBQueryStruct


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

    /** Set query function information for subfield.
     *
     * The default is to use the database passed in the call to openDB().
     * Passing in a spatial database in this function overrised use of the
     * default database.
     *
     * @param[in] subfield Name of subfield.
     * @param[in] fn Query function to use for subfield.
     * @param[in] db Spatial database to query.
     */
    void queryFn(const char* subfield,
                 const queryfn_type fn,
                 spatialdata::spatialdb::SpatialDB* db =NULL);

    /** Get query function for subfield.
     *
     * @param subfield Name of subfield.
     * @return Query function used for subfield.
     */
    const queryfn_type queryFn(const char* subfield) const;

    /** Get spatial database used to get values for subfield.
     *
     * @param subfield Name of subfield.
     * @return Spatial database to use in query (NULL if using default database passed in openDB()).
     */
    const spatialdata::spatialdb::SpatialDB* queryDB(const char* subfield) const;

    /// Get array of query functions.
    queryfn_type* functions(void) const;

    /// Get array of pointers to contexts.
    const DBQueryContext* const* contextPtrs(void) const;

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

    /** Close spatial database query for setting values in field.
     *
     * @param db Spatial database to query.
     */
    void closeDB(spatialdata::spatialdb::SpatialDB* db);


    /** Generic query of values from spatial database.
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
    PetscErrorCode dbQueryGeneric(PylithInt dim,
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


    // PRIVATE TYPEDEFS /////////////////////////////////////////////////////
private:

    typedef std::map<std::string, queryfn_type> queryfn_map_type;
    typedef std::map<std::string, spatialdata::spatialdb::SpatialDB*> querydb_map_type;


    // PRIVATE MEMBERS //////////////////////////////////////////////////////
private:

    const pylith::topology::Field& _field; ///< Field associated with query.

    queryfn_map_type _queryFns;
    querydb_map_type _queryDBs;

    queryfn_type* _functions; ///< Functions implementing queries.
    DBQueryContext* _contexts; ///< Contexts for performing query for each subfield.
    DBQueryContext** _contextPtrs; ///< Array of pointers to contexts.

    // NOT IMPLEMENTED //////////////////////////////////////////////////////
private:

    FieldQuery(const FieldQuery&); ///< Not implemented
    const FieldQuery& operator=(const FieldQuery&); ///< Not implemented

}; // FieldQuery

#endif // pylith_topology_fieldquery_hh


// End of file
