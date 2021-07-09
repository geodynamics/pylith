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

/** @file libsrc/materials/Quuery.hh
 *
 */

#if !defined(pylith_materials_query_hh)
#define pylith_materials_query_hh

#include "pylith/materials/materialsfwd.hh" // forward declarations

#include "pylith/feassemble/feassemblefwd.hh" // USES AuxiliaryFactory
#include "spatialdata/spatialdb/spatialdbfwd.hh" // USES SpatialDB

#include <cstddef> // USES size_t

class pylith::materials::Query {
    friend class TestQuery; // unit testing

    // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    /** Setup subfield query in auxiliary factory for shear modulus from density and Vs.
     *
     * @param[in] subfieldName Name for shear modulus subfield.
     * @param[inout] factory Auxiliary factory associated with shear modulus subfield.
     */
    static
    void shearModulusFromVM(const char* subfieldName,
                            pylith::feassemble::AuxiliaryFactory* factory);

    /** Setup subfield query in auxiliary factory for bulk modulus from density, Vs, and Vp.
     *
     * @param[in] subfieldName Name for bulk modulus subfield.
     * @param[inout] factory Auxiliary factory associated with shear modulus subfield.
     */
    static
    void bulkModulusFromVM(const char* subfieldName,
                           pylith::feassemble::AuxiliaryFactory* factory);

    /** Setup subfield query in auxiliary factory for Maxwell time from density, Vs, and viscosity.
     *
     * @param[in] subfieldName Name for Maxwell time subfield.
     * @param[inout] factory Auxiliary factory associated with shear modulus subfield.
     */
    static
    void maxwellTimeFromVM(const char* subfieldName,
                           pylith::feassemble::AuxiliaryFactory* factory);

    /** Setup subfield query in auxiliary factory for generalized Maxwell times from density, Vs, and viscosities.
     *
     * @param[in] subfieldName Name for generalized Maxwelltime subfield.
     * @param[inout] factory Auxiliary factory associated with shear modulus subfield.
     */
    static
    void generalizedMaxwellTimesFromVM(const char* subfieldName,
                                       pylith::feassemble::AuxiliaryFactory* factory);

    /** Setup subfield query in auxiliary factory for generalized Maxwell shear modulus ratios.
     *
     * @param[in] subfieldName Name for generalized Maxwelltime subfield.
     * @param[inout] factory Auxiliary factory associated with shear modulus subfield.
     */
    static
    void generalizedMaxwellShearModulusRatiosFromVM(const char* subfieldName,
                                                    pylith::feassemble::AuxiliaryFactory* factory);

    /** Setup subfield query in auxiliary factory for gravity field from GravityField spatial database.
     *
     * @param[in] subfieldName Name for shear modulus subfield.
     * @param[inout] factory Auxiliary factory associated with shear modulus subfield.
     * @param[in] gravityField Spatial database for gravity field.
     * @param[in] Spatial dimension of problem.
     */
    static
    void gravityFieldFromDB(const char* subfieldName,
                            pylith::feassemble::AuxiliaryFactory* factory,
                            spatialdata::spatialdb::GravityField* gravityField,
                            const size_t spaceDim);

    /** Setup subfield query in auxiliary factory for bulk modulus from  gravity field from input parameters.
     *
     * @param[in] subfieldName Name for shear modulus subfield.
     * @param[inout] factory Auxiliary factory associated with shear modulus subfield.
     */
    static
    void biotModulusFromInput(const char* subfieldName,
                            pylith::feassemble::AuxiliaryFactory* factory);


    // NOT IMPLEMENTED //////////////////////////////////////////////////////
private:

    Query(void); ///< Not implemented
    ~Query(void); ///< Not implemented
    Query(const Query&); ///< Not implemented
    const Query& operator=(const Query&); ///< Not implemented

}; // Query

#endif // pylith_materials_query_hh

// End of file
