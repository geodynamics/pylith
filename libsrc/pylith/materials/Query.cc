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

#include <portinfo>

#include "pylith/materials/Query.hh"

#include "pylith/topology/FieldQuery.hh" // USES DBQueryContext
#include "pylith/feassemble/AuxiliaryFactory.hh" // USES AuxiliaryFactory
#include "pylith/utils/error.hh" // USES PYLITH_METHOD_BEGIN/END, PYLITH_ERROR_RETURN
#include "pylith/utils/array.hh" // USES scalar_array, int_array
#include "pylith/utils/types.hh" // USES PylithScalar
#include "pylith/utils/constdefs.h" // USES PYLITH_MAXSCALAR

#include "spatialdata/spatialdb/GravityField.hh" // USES GravityField

namespace pylith {
    namespace materials {
        class _Query {
            /** Conversion functions to convert from spatial database values to auxiliary subfield components.
             *
             * Order of spatial database values must match the order given in the Query factory functions.
             *
             * Interface
             *
             * @param[out] valueSubfield Value for subfield.
             * @param[in] numComponents Number of components for subfield value.
             * @param[in] dbValues Array of values from spatial database query.
             * @param[in] dbIndices Indices of values from spatial database to use for computing subfield values.
             * @returns PETSc error code (0 for success).
             */
public:

            static
            std::string vmToShearModulus(PylithScalar valueSubfield[],
                                         const PylithInt numComponents,
                                         const pylith::scalar_array dbValues,
                                         const pylith::int_array dbIndices);

            static
            std::string vmToBulkModulus(PylithScalar valueSubfield[],
                                        const PylithInt numComponents,
                                        const pylith::scalar_array dbValues,
                                        const pylith::int_array dbIndices);

            static
            std::string vmToMaxwellTime(PylithScalar valueSubfield[],
                                        const PylithInt numComponents,
                                        const pylith::scalar_array dbValues,
                                        const pylith::int_array dbIndices);

            static
            std::string vmToGeneralizedMaxwellTimes(PylithScalar valueSubfield[],
                                                    const PylithInt numComponents,
                                                    const pylith::scalar_array dbValues,
                                                    const pylith::int_array dbIndices);

            static
            std::string vmToGeneralizedMaxwellShearModulusRatios(PylithScalar valueSubfield[],
                                                                 const PylithInt numComponents,
                                                                 const pylith::scalar_array dbValues,
                                                                 const pylith::int_array dbIndices);

            static
            std::string dbToGravityField(PylithScalar valueSubfield[],
                                         const PylithInt numComponents,
                                         const pylith::scalar_array dbValues,
                                         const pylith::int_array dbIndices);

            static
            std::string inputToBiotModulus(PylithScalar valueSubfield[],
                                           const PylithInt numComponents,
                                           const pylith::scalar_array dbValues,
                                           const pylith::int_array dbIndices);

        }; // _Query
    } // materials
} // pylith

// ----------------------------------------------------------------------
// Setup subfield query in auxiliary factory for shear modulus from density and Vs.
void
pylith::materials::Query::shearModulusFromVM(const char* subfieldName,
                                             pylith::feassemble::AuxiliaryFactory* factory) {
    const size_t numDBValues = 2;
    const char* dbValues[numDBValues] = { "density", "vs" };

    assert(factory);
    factory->setSubfieldQuery(subfieldName, dbValues, numDBValues, _Query::vmToShearModulus);
} // shearModulusFromVM


// ----------------------------------------------------------------------
// Setup subfield query in auxiliary factory for bulk modulus from density, Vs, and Vp.
void
pylith::materials::Query::bulkModulusFromVM(const char* subfieldName,
                                            pylith::feassemble::AuxiliaryFactory* factory) {
    const size_t numDBValues = 3;
    const char* dbValues[numDBValues] = { "density", "vs", "vp" };

    assert(factory);
    factory->setSubfieldQuery(subfieldName, dbValues, numDBValues, _Query::vmToBulkModulus);
} // bulkModulusFromVM


// ----------------------------------------------------------------------
// Setup subfield query in auxiliary factory for Maxwell time from density, Vs, and viscosity.
void
pylith::materials::Query::maxwellTimeFromVM(const char* subfieldName,
                                            pylith::feassemble::AuxiliaryFactory* factory) {
    const size_t numDBValues = 3;
    const char* dbValues[numDBValues] = { "density", "vs", "viscosity" };

    assert(factory);
    factory->setSubfieldQuery(subfieldName, dbValues, numDBValues, _Query::vmToMaxwellTime);
} // maxwellTimeFromVM


// ----------------------------------------------------------------------
// Setup subfield query in auxiliary factory for generalized Maxwell times from density, Vs, and viscosities.
void
pylith::materials::Query::generalizedMaxwellTimesFromVM(const char* subfieldName,
                                                        pylith::feassemble::AuxiliaryFactory* factory) {
    const size_t numDBValues = 8;
    const char* dbValues[numDBValues] = {
        "density",
        "vs",
        "viscosity_1",
        "viscosity_2",
        "viscosity_3",
        "shear_modulus_ratio_1",
        "shear_modulus_ratio_2",
        "shear_modulus_ratio_3",
    };

    assert(factory);
    factory->setSubfieldQuery(subfieldName, dbValues, numDBValues, _Query::vmToGeneralizedMaxwellTimes);
} // generalizedMaxwellTimesFromVM


// ----------------------------------------------------------------------
// Setup subfield query in auxiliary factory for generalized Maxwell shear modulus ratios.
void
pylith::materials::Query::generalizedMaxwellShearModulusRatiosFromVM(const char* subfieldName,
                                                                     pylith::feassemble::AuxiliaryFactory* factory) {
    const size_t numDBValues = 3;
    const char* dbValues[numDBValues] = {
        "shear_modulus_ratio_1",
        "shear_modulus_ratio_2",
        "shear_modulus_ratio_3",
    };

    assert(factory);
    factory->setSubfieldQuery(subfieldName, dbValues, numDBValues, _Query::vmToGeneralizedMaxwellShearModulusRatios);
} // generalizedMaxwellShearModulusRatiosFromVM


// ----------------------------------------------------------------------
// Setup subfield query in auxiliary factory for gravity field from GravityField spatial database.
void
pylith::materials::Query::gravityFieldFromDB(const char* subfieldName,
                                             pylith::feassemble::AuxiliaryFactory* factory,
                                             spatialdata::spatialdb::GravityField* gravityField,
                                             const size_t spaceDim) {
    const size_t numDBValues = 3;
    const char* dbValues[numDBValues] = {
        "gravity_field_x",
        "gravity_field_y",
        "gravity_field_z"
    };

    assert(factory);
    factory->setSubfieldQuery(subfieldName, dbValues, spaceDim, _Query::dbToGravityField, gravityField);
} // gravityFieldFromDB

// ----------------------------------------------------------------------
// Setup subfield query in auxiliary factory for bulk modulus from density, Vs, and Vp.
void
pylith::materials::Query::biotModulusFromInput(const char* subfieldName,
                                            pylith::feassemble::AuxiliaryFactory* factory) {
    const size_t numDBValues = 4;
    const char* dbValues[numDBValues] = {
        "fluid_bulk_modulus",
        "solid_bulk_modulus",
        "biot_coefficient",
        "porosity"
    };

    assert(factory);
    factory->setSubfieldQuery(subfieldName, dbValues, numDBValues, _Query::inputToBiotModulus);
} // biotModulusFromInput


// ----------------------------------------------------------------------
// Compute shear modules from density and Vs.
std::string
pylith::materials::_Query::vmToShearModulus(PylithScalar valueSubfield[],
                                            const PylithInt numComponents,
                                            const pylith::scalar_array dbValues,
                                            const pylith::int_array dbIndices) {
    PYLITH_METHOD_BEGIN;

    const size_t _numComponents = 1;

    assert(valueSubfield);
    assert(_numComponents == size_t(numComponents));
    assert(2 == dbIndices.size());

    const size_t i_density = 0;assert(dbIndices[i_density] < dbValues.size());
    const size_t i_vs = 1;assert(dbIndices[i_vs] < dbValues.size());
    const PylithScalar density = dbValues[dbIndices[i_density]];
    const PylithScalar vs = dbValues[dbIndices[i_vs]];
    valueSubfield[0] = density * vs * vs;

    bool valuesOkay = true;
    std::ostringstream msg;
    if (density <= 0) {
        valuesOkay = false;
        msg << "Found negative density (" << density << ").";
    } // if
    if (vs <= 0) {
        valuesOkay = false;
        msg << "Found negative shear wave speed (" << vs << ").";
    } // if

    PYLITH_METHOD_RETURN(msg.str());
} // vmToShearModulus


// ----------------------------------------------------------------------
// Compute bulk modules from density, Vs, and Vp.
std::string
pylith::materials::_Query::vmToBulkModulus(PylithScalar valueSubfield[],
                                           const PylithInt numComponents,
                                           const pylith::scalar_array dbValues,
                                           const pylith::int_array dbIndices) {
    PYLITH_METHOD_BEGIN;

    const size_t _numComponents = 1;

    assert(valueSubfield);
    assert(_numComponents == size_t(numComponents));
    assert(3 == dbIndices.size());

    const size_t i_density = 0;assert(dbIndices[i_density] < dbValues.size());
    const size_t i_vs = 1;assert(dbIndices[i_vs] < dbValues.size());
    const size_t i_vp = 2;assert(dbIndices[i_vp] < dbValues.size());
    const PylithScalar density = dbValues[dbIndices[i_density]];
    const PylithScalar vs = dbValues[dbIndices[i_vs]];
    const PylithScalar vp = dbValues[dbIndices[i_vp]];
    valueSubfield[0] = density * (vp*vp - 4.0/3.0*vs*vs);

    bool valuesOkay = true;
    std::ostringstream msg;
    if (density <= 0) {
        valuesOkay = false;
        msg << "Found nonpositive density (" << density << ").";
    } // if
    if (vs < 0) {
        valuesOkay = false;
        msg << "Found negative shear wave speed (" << vs << ").";
    } // if
    if (vp <= 0) {
        valuesOkay = false;
        msg << "Found nonpositive dilatational wave speed (" << vp << ").";
    } // if

    PYLITH_METHOD_RETURN(msg.str());
} // vmToBulkModulus


// ----------------------------------------------------------------------
// Compute Maxwell time from from density, Vs, and viscosity.
std::string
pylith::materials::_Query::vmToMaxwellTime(PylithScalar valueSubfield[],
                                           const PylithInt numComponents,
                                           const pylith::scalar_array dbValues,
                                           const pylith::int_array dbIndices) {
    PYLITH_METHOD_BEGIN;

    const size_t _numComponents = 1;

    assert(valueSubfield);
    assert(_numComponents == size_t(numComponents));
    assert(3 == dbIndices.size());

    const size_t i_density = 0;assert(dbIndices[i_density] < dbValues.size());
    const size_t i_vs = 1;assert(dbIndices[i_vs] < dbValues.size());
    const size_t i_viscosity = 2;assert(dbIndices[i_viscosity] < dbValues.size());
    const PylithScalar density = dbValues[dbIndices[i_density]];
    const PylithScalar vs = dbValues[dbIndices[i_vs]];
    const PylithScalar viscosity = dbValues[dbIndices[i_viscosity]];
    const PylithScalar shearModulus = density * vs * vs;
    valueSubfield[0] = viscosity / shearModulus;

    bool valuesOkay = true;
    std::ostringstream msg;
    if (density <= 0) {
        valuesOkay = false;
        msg << "Found negative density (" << density << ").";
    } // if
    if (vs <= 0) {
        valuesOkay = false;
        msg << "Found negative shear wave speed (" << vs << ").";
    } // if
    if (viscosity <= 0) {
        valuesOkay = false;
        msg << "Found nonpositive viscosity (" << viscosity << ").";
    } // if

    PYLITH_METHOD_RETURN(msg.str());
} // vmToMaxwellTime


// ----------------------------------------------------------------------
// Compute Maxwell time for generalized Maxwell model (3 elements).
std::string
pylith::materials::_Query::vmToGeneralizedMaxwellTimes(PylithScalar valueSubfield[],
                                                       const PylithInt numComponents,
                                                       const pylith::scalar_array dbValues,
                                                       const pylith::int_array dbIndices) {
    PYLITH_METHOD_BEGIN;

    const size_t _numComponents = 3;

    assert(valueSubfield);
    assert(_numComponents == size_t(numComponents));
    assert(8 == dbIndices.size());

    const size_t i_density = 0;assert(dbIndices[i_density] < dbValues.size());
    const size_t i_vs = 1;assert(dbIndices[i_vs] < dbValues.size());
    const size_t i_viscosity = 2;assert(dbIndices[i_viscosity]+2 < dbValues.size());
    const size_t i_shearModulusRatio = 5;assert(dbIndices[i_shearModulusRatio]+2 < dbValues.size());

    const PylithScalar density = dbValues[dbIndices[i_density]];
    const PylithScalar vs = dbValues[dbIndices[i_vs]];
    const PylithScalar viscosity1 = dbValues[dbIndices[i_viscosity+0]];
    const PylithScalar viscosity2 = dbValues[dbIndices[i_viscosity+1]];
    const PylithScalar viscosity3 = dbValues[dbIndices[i_viscosity+2]];
    const PylithScalar shearModulusRatio1 = dbValues[dbIndices[i_shearModulusRatio+0]];
    const PylithScalar shearModulusRatio2 = dbValues[dbIndices[i_shearModulusRatio+1]];
    const PylithScalar shearModulusRatio3 = dbValues[dbIndices[i_shearModulusRatio+2]];

    const PylithScalar shearModulus = density * vs * vs;

    const PylithScalar shearModulus1 = shearModulusRatio1 * shearModulus;
    valueSubfield[0] = (shearModulus1 > 0.0) ? viscosity1 / shearModulus1 : PYLITH_MAXSCALAR;

    const PylithReal shearModulus2 = shearModulusRatio2 * shearModulus;
    valueSubfield[1] = (shearModulus2 > 0.0) ? viscosity2 / shearModulus2 : PYLITH_MAXSCALAR;

    const PylithReal shearModulus3 = shearModulusRatio3 * shearModulus;
    valueSubfield[2] = (shearModulus3 > 0.0) ? viscosity3 / shearModulus3 : PYLITH_MAXSCALAR;

    bool valuesOkay = true;
    std::ostringstream msg;
    if (density <= 0) {
        valuesOkay = false;
        msg << "Found negative density (" << density << ").";
    } // if
    if (vs <= 0) {
        valuesOkay = false;
        msg << "Found negative shear wave speed (" << vs << ").";
    } // if
    if (viscosity1 <= 0) {
        valuesOkay = false;
        msg << "Found nonpositive viscosity 1 (" << viscosity1 << ").";
    } // if
    if (viscosity2 <= 0) {
        valuesOkay = false;
        msg << "Found nonpositive viscosity 2 (" << viscosity2 << ").";
    } // if
    if (viscosity3 <= 0) {
        valuesOkay = false;
        msg << "Found nonpositive viscosity 3 (" << viscosity3 << ").";
    } // if
    if (shearModulusRatio1 <= 0) {
        valuesOkay = false;
        msg << "Found negative shear modulus ratio 1 (" << shearModulusRatio1 << ").";
    } // if
    if (shearModulusRatio2 <= 0) {
        valuesOkay = false;
        msg << "Found negative shear modulus ratio 2 (" << shearModulusRatio2 << ").";
    } // if
    if (shearModulusRatio3 <= 0) {
        valuesOkay = false;
        msg << "Found negative shear modulus ratio 3 (" << shearModulusRatio3 << ").";
    } // if

    const double ratioSum = shearModulusRatio1 + shearModulusRatio2 + shearModulusRatio3;
    if (ratioSum > 1) {
        valuesOkay = false;
        msg << "Shear ratio sum greater than one (" << ratioSum << ").";
        msg << " Shear ratios are " << shearModulusRatio1 << ", " << shearModulusRatio2 << ", " << shearModulusRatio3
            << ".";
    } // if

    PYLITH_METHOD_RETURN(msg.str());
} // vmToGeneralizedMaxwellTimes


// ----------------------------------------------------------------------
// Compute Maxwell time for generalized Maxwell model (3 elements).
std::string
pylith::materials::_Query::vmToGeneralizedMaxwellShearModulusRatios(PylithScalar valueSubfield[],
                                                                    const PylithInt numComponents,
                                                                    const pylith::scalar_array dbValues,
                                                                    const pylith::int_array dbIndices) {
    PYLITH_METHOD_BEGIN;

    const size_t _numComponents = 3;

    assert(valueSubfield);
    assert(_numComponents == size_t(numComponents));
    assert(3 == dbIndices.size());

    const size_t i_shearModulusRatio = 0;assert(dbIndices[i_shearModulusRatio+2] < dbValues.size());

    const PylithScalar shearModulusRatio1 = valueSubfield[0] = dbValues[dbIndices[i_shearModulusRatio+0]];
    const PylithScalar shearModulusRatio2 = valueSubfield[1] = dbValues[dbIndices[i_shearModulusRatio+1]];
    const PylithScalar shearModulusRatio3 = valueSubfield[2] = dbValues[dbIndices[i_shearModulusRatio+2]];

    bool valuesOkay = true;
    std::ostringstream msg;
    if (shearModulusRatio1 <= 0) {
        valuesOkay = false;
        msg << "Found negative shear modulus ratio 1 (" << shearModulusRatio1 << ").";
    } // if
    if (shearModulusRatio2 <= 0) {
        valuesOkay = false;
        msg << "Found negative shear modulus ratio 2 (" << shearModulusRatio2 << ").";
    } // if
    if (shearModulusRatio3 <= 0) {
        valuesOkay = false;
        msg << "Found negative shear modulus ratio 3 (" << shearModulusRatio3 << ").";
    } // if
    const double ratioSum = shearModulusRatio1 + shearModulusRatio2 + shearModulusRatio3;
    if (ratioSum > 1) {
        valuesOkay = false;
        msg << "Shear ratio sum greater than one (" << ratioSum << ").";
        msg << " Shear ratios are " << shearModulusRatio1 << ", " << shearModulusRatio2 << ", " << shearModulusRatio3
            << ".";
    } // if

    PYLITH_METHOD_RETURN(msg.str());
} // vmToGeneralizedMaxwellShearModulusRatios


// ----------------------------------------------------------------------
// Compute Maxwell time for generalized Maxwell model (3 elements).
std::string
pylith::materials::_Query::dbToGravityField(PylithScalar valueSubfield[],
                                            const PylithInt numComponents,
                                            const pylith::scalar_array dbValues,
                                            const pylith::int_array dbIndices) {
    PYLITH_METHOD_BEGIN;

    const size_t spaceDim = dbIndices.size();

    assert(valueSubfield);
    assert(spaceDim == size_t(numComponents));

    for (size_t i = 0; i < spaceDim; ++i) {
        assert(dbIndices[i] < dbValues.size());
        valueSubfield[i] = dbValues[dbIndices[i]];
    } // for

    PylithScalar mag = 0.0;
    for (size_t i = 0; i < spaceDim; ++i) {
        mag += valueSubfield[i] * valueSubfield[i];
    } // for
    const PylithReal tolerance = 1.0e-6;
    std::ostringstream msg;
    if (mag < tolerance) {
        msg << "Found near zero magnitude (" << mag << ") for gravity field vector (";
        for (size_t i = 0; i < spaceDim; ++i) {
            msg << "  " << valueSubfield[i];
        } // for
        msg << ".";
    } // if

    PYLITH_METHOD_RETURN(msg.str());
} // dbToGravityField

// ----------------------------------------------------------------------
// Compute Biot's modulus from Biot's coefficient, solid grain bulk moduls,
// fluid bulk modulus, and porosity
std::string
pylith::materials::_Query::inputToBiotModulus(PylithScalar valueSubfield[],
                                            const PylithInt numComponents,
                                            const pylith::scalar_array dbValues,
                                            const pylith::int_array dbIndices) {
    PYLITH_METHOD_BEGIN;

    const size_t _numComponents = 1;

    assert(valueSubfield);
    assert(_numComponents == size_t(numComponents));
    assert(4 == dbIndices.size());

    const size_t i_fluid_bulk_modulus = 0;assert(dbIndices[i_fluid_bulk_modulus] < dbValues.size());
    const size_t i_solid_bulk_modulus = 1;assert(dbIndices[i_solid_bulk_modulus] < dbValues.size());
    const size_t i_biot_coefficient = 2;assert(dbIndices[i_biot_coefficient] < dbValues.size());
    const size_t i_porosity = 3;assert(dbIndices[i_porosity] < dbValues.size());

    const PylithScalar fluid_bulk_modulus = dbValues[dbIndices[i_fluid_bulk_modulus]];
    const PylithScalar solid_bulk_modulus = dbValues[dbIndices[i_solid_bulk_modulus]];
    const PylithScalar biot_coefficient = dbValues[dbIndices[i_biot_coefficient]];
    const PylithScalar porosity = dbValues[dbIndices[i_porosity]];

    valueSubfield[0] = 1.0 / ( porosity / fluid_bulk_modulus + (biot_coefficient - porosity) / solid_bulk_modulus );

    bool valuesOkay = true;
    std::ostringstream msg;
    if (porosity < 0) {
        valuesOkay = false;
        msg << "Found negative porosity (" << porosity << ").";
    } // if
    if (biot_coefficient <= 0) {
        valuesOkay = false;
        msg << "Found negative biot coefficient (" << biot_coefficient << ").";
    } // if

    // Debug
    PylithScalar biot_modulus = 1.0 / ( porosity / fluid_bulk_modulus + (biot_coefficient - porosity) / solid_bulk_modulus );
    if (biot_modulus <= 0) {
        valuesOkay = false;
        msg << "biot modulus (" << biot_modulus << ") wrong. Kfl: " << fluid_bulk_modulus << " Ksg: " << solid_bulk_modulus << " phi: " << porosity << " alpha: " << biot_coefficient;
    } // if

    PYLITH_METHOD_RETURN(msg.str());
} // inputToBiotModulus

// End of file
