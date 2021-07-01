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

/** @file libsrc/bc/AuxiliaryFactory.hh
 *
 * @brief C++ helper class for setting up auxiliary subfields.
 */

#if !defined(pylith_feassemble_auxiliaryfactory_hh)
#define pylith_feassemble_auxiliaryfactory_hh

#include "feassemblefwd.hh" // forward declarations
#include "pylith/topology/FieldFactory.hh" // ISA FieldFactory

#include "pylith/topology/FieldBase.hh" // USES FieldBase::Discretization
#include "pylith/topology/FieldQuery.hh" // USES FieldQuery::queryfn_type

#include "spatialdata/spatialdb/spatialdbfwd.hh" // USES SpatialDB
#include "spatialdata/units/unitsfwd.hh" // HOLDSA Normalizer

class pylith::feassemble::AuxiliaryFactory : public pylith::topology::FieldFactory {
    friend class TestAuxiliaryFactory; // unit testing

    // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    /// Default constructor.
    AuxiliaryFactory(void);

    /// Destructor.
    virtual ~AuxiliaryFactory(void);

    /** Set spatial database for filling auxiliary subfields.
     *
     * @param[in] value Pointer to database.
     */
    void setQueryDB(spatialdata::spatialdb::SpatialDB* value);

    /** Get spatial database for filling auxiliary subfields.
     *
     * @returns Pointer to database.
     */
    const spatialdata::spatialdb::SpatialDB* getQueryDB(void) const;

    /** Initialize factory for setting up auxiliary subfields.
     *
     * @param[inout] field Auxiliary field for which subfields are to be created.
     * @param[in] normalizer Scales for nondimensionalization.
     * @param[in] spaceDim Spatial dimension of problem.
     * @param[in] defaultDescription Default description for new subfields.
     */
    void initialize(pylith::topology::Field* field,
                    const spatialdata::units::Nondimensional& normalizer,
                    const int spaceDim,
                    const pylith::topology::FieldBase::Description* defaultDescription=NULL);

    /// Set subfield values using spatial database.
    void setValuesFromDB(void);

    /** Set query function for subfield.
     *
     * @param[in] subfieldName Name of subfield.
     * @param[in] namesDBValues Array of names of values to use from spatial database.
     * @param[in] numDBValues Size of names array.
     * @param[in] convertFn Function to convert spatial database values to subfield values.
     * @param[in] db Spatial database to query.
     */
    void setSubfieldQuery(const char* subfieldName,
                          const char* namesDBValues[]=NULL,
                          const size_t numDBValues=0,
                          pylith::topology::FieldQuery::convertfn_type convertFn=NULL,
                          spatialdata::spatialdb::SpatialDB* db=NULL);

    // PROTECTED MEMBERS ///////////////////////////////////////////////////////////////////////////////////////////////
protected:

    /** Database of values for filling subfields.
     * Auxiliary subfields are not filled if NULL.
     *
     * Currently, this is the only way to populate the auxiliary subfields.
     */
    spatialdata::spatialdb::SpatialDB* _queryDB;

    /// Field query for filling subfield values via spatial database.
    pylith::topology::FieldQuery* _fieldQuery;

    // NOT IMPLEMENTED /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    AuxiliaryFactory(const AuxiliaryFactory &); ///< Not implemented.
    const AuxiliaryFactory& operator=(const AuxiliaryFactory&); ///< Not implemented

}; // class AuxiliaryFactory

#endif // pylith_feassemble_auxiliaryfactory_hh

// End of file
