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
 * @file libsrc/testing/FieldTester.hh
 *
 * @brief General routines for C++ unit tests related to Field.
 */

#if !defined(pylith_testing_fieldtester_hh)
#define pylith_testing_fieldtester_hh

#include "pylith/testing/testingfwd.hh" // forward declarations

#include "pylith/topology/Field.hh" // USES Field::SubfieldInfo
#include "pylith/utils/petscfwd.h" // USES PetscFE

#include "spatialdata/spatialdb/spatialdbfwd.hh" // USES SpatialDB

class pylith::testing::FieldTester {
    // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    /** Check to make sure field matches spatial database.
     *
     * @param[in] field Field to check.
     * @param[in] fieldDB Spatial database describing field.
     * @param[in] lengthScale Length scale for nondimensionalization.
     * @returns L2 norm of difference between field and spatial database.
     */
    static
    PylithReal checkFieldWithDB(const pylith::topology::Field& field,
                                spatialdata::spatialdb::SpatialDB* fieldDB,
                                const PylithReal lengthScale);

    /** Check that subfield information in field test subject matches expected subfield.
     *
     * @param field Field with subfields created by factory.
     * @param infoE Expected subfield info.
     */
    static
    void checkSubfieldInfo(const pylith::topology::Field& field,
                           const pylith::topology::Field::SubfieldInfo& infoE);

    // NOT IMPLEMENTED //////////////////////////////////////////////////////
private:

    FieldTester(void); ///< Not implemented.
    FieldTester(const FieldTester&); ///< Not implemented.
    const FieldTester& operator=(const FieldTester&); ///< Not implemented.

}; // FieldTester

#endif // pylith_testing_fieldtester_hh

// End of file
