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

/** @file libsrc/bc/AuxiliaryFactory.hh
 *
 * @brief C++ helper class for setting up auxiliary subfields.
 */

#if !defined(pylith_feassemble_auxiliaryfactory_hh)
#define pylith_feassemble_auxiliaryfactory_hh

#include "feassemblefwd.hh" // forward declarations
#include "pylith/utils/GenericComponent.hh" // ISA GenericComponent

#include "pylith/topology/FieldBase.hh" // USES FieldBase::Discretization
#include "pylith/topology/FieldQuery.hh" // USES FieldQuery::queryfn_type

#include "spatialdata/spatialdb/spatialdbfwd.hh" // USES SpatialDB
#include "spatialdata/units/unitsfwd.hh" // HOLDSA Normalizer

class pylith::feassemble::AuxiliaryFactory : public pylith::utils::GenericComponent {
    friend class TestAuxiliaryFactory; // unit testing

    // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    /// Default constructor.
    AuxiliaryFactory(void);

    /// Destructor.
    ~AuxiliaryFactory(void);

    /** Set spatial database for filling auxiliary subfields.
     *
     * @param[in] value Pointer to database.
     */
    void setQueryDB(spatialdata::spatialdb::SpatialDB* value);

    /** Get spatial database for filling auxiliary subfields.
     *
     * @returns Pointer to database.
     */
    const spatialdata::spatialdb::SpatialDB* getQueryDB(void);

    /** Set discretization information for auxiliary subfield.
     *
     * @param[in] subfieldName Name of auxiliary subfield.
     * @param[in] basisOrder Polynomial order for basis.
     * @param[in] quadOrder Order of quadrature rule.
     * @param[in] dimension Dimension of points for discretization.
     * @param[in] isBasisContinuous True if basis is continuous.
     * @param[in] feSpace Finite-element space.
     */
    void setSubfieldDiscretization(const char* subfieldName,
                                   const int basisOrder,
                                   const int quadOrder,
                                   const int dimension,
                                   const bool isBasisContinuous,
                                   const pylith::topology::FieldBase::SpaceEnum feSpace);

    /** Get discretization information for subfield.
     *
     * @param[in] subfieldName Name of subfield.
     * @return Discretization information for auxiliary subfield. If
     * discretization information was not set, then use "default".
     */
    const pylith::topology::FieldBase::Discretization& getSubfieldDiscretization(const char* subfieldName) const;

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

    // PROTECTED METHODS ///////////////////////////////////////////////////////////////////////////////////////////////
protected:

    /** Set query function for subfield.
     *
     * @param[in] subfieldName Name of subfield.
     * @param[in] fn Function used for query spatial database.
     * @param[in] db Spatial database to query.
     */
    void _setSubfieldQueryFn(const char* subfieldName,
                             pylith::topology::FieldQuery::queryfn_type,
                             spatialdata::spatialdb::SpatialDB* db = NULL);

    // PROTECTED MEMBERS ///////////////////////////////////////////////////////////////////////////////////////////////
protected:

    pylith::topology::Field* _field; ///< Auxiliary field.
    pylith::topology::FieldBase::discretizations_map _subfieldDiscretizations; ///< Discretization for each subfield.
    pylith::topology::FieldBase::Description* _defaultDescription; ///< Description for default subfield.
    spatialdata::units::Nondimensional* _normalizer; ///< Scales for nondimensionalization.
    int _spaceDim; ///< Spatial dimension.

    /** Database of values for filling subfields.
     * Auxiliary subfields are not filled if NULL.
     *
     * Currently, this is the only way to fill subfields.
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
