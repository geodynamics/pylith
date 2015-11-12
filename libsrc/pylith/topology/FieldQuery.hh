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
public :

  /// Function prototype for query functions.
  typedef PetscErrorCode (*queryfn_type)(PetscInt, PetscReal, const PetscReal[], PetscInt, PetscScalar*, void*);

// PUBLIC STRUCT ////////////////////////////////////////////////////////
public :

  /// Context for spatial database queries.
  struct DBQueryContext {
    spatialdata::spatialdb::SpatialDB* db; ///< Spatial database.
    const spatialdata::geocoords::CoordSys* cs; ///< Coordinate system of point locations.
    PylithReal lengthScale; ///< Length scale for dimensionalizing coordinates.
    PylithReal valueScale; ///< Scale for dimensionalizing values for subfield.
  }; // dbQueryStruct
  

// PUBLIC MEMBERS ///////////////////////////////////////////////////////
public :

  /** Default constructor.
   *
   * @param field Field associated with query.
   */
  FieldQuery(const Field& field);

  /// Destructor.
  ~FieldQuery(void);

  /// Deallocate PETSc and local data structures.
  void deallocate(void);

  /** Set query function for subfield.
   *
   * @param subfield Name of subfield.
   * @param fn Query function to use for subfield.
   */
  void queryFn(const char* subfield,
	       const queryfn_type fn);

  /** Get query function for subfield.
   *
   * @param subfield Name of subfield.
   * @return Query function used for subfield.
   */
  const queryfn_type queryFn(const char* subfield) const;

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
  

// PRIVATE TYPEDEFS /////////////////////////////////////////////////////
private :

  typedef std::map<std::string, queryfn_type> queryfn_map_type;


// PRIVATE MEMBERS //////////////////////////////////////////////////////
private :

  const pylith::topology::Field& _field; ///< Field associated with query.

  queryfn_map_type _queryFns;

  queryfn_type* _functions; ///< Functions implementing queries.
  DBQueryContext* _contexts; ///< Contexts for performing query for each subfield.
  DBQueryContext** _contextPtrs; ///< Array of pointers to contexts.

// NOT IMPLEMENTED //////////////////////////////////////////////////////
private :

  FieldQuery(const FieldQuery&); ///< Not implemented
  const FieldQuery& operator=(const FieldQuery&); ///< Not implemented

}; // FieldQuery

#endif // pylith_topology_fieldquery_hh


// End of file 
