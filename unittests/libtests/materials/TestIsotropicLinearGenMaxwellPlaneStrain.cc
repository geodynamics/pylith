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
// Copyright (c) 2010-2015 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "TestIsotropicLinearGenMaxwellPlaneStrain.hh" // Implementation of class methods

#include "pylith/materials/IsotropicLinearGenMaxwellPlaneStrain.hh" // USES IsotropicLinearGenMaxwellPlaneStrain
#include "pylith/materials/Query.hh" // USES Query

#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/Fields.hh" // USES Fields
#include "pylith/topology/FieldQuery.hh" // USES FieldQuery
#include "pylith/problems/SolutionFactory.hh" // USES SolutionFactory

#include "pylith/utils/error.hh" // USES PYLITH_METHOD_BEGIN/END
#include "pylith/utils/journals.hh" // USES JournalingComponent

#include "spatialdata/spatialdb/UserFunctionDB.hh" // USES UserFunctionDB
#include "spatialdata/geocoords/CSCart.hh" // USES CSCart
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional


// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::materials::TestIsotropicLinearGenMaxwellPlaneStrain::setUp(void) {
    TestMaterial::setUp();
    _mymaterial = new IsotropicLinearGenMaxwellPlaneStrain(); CPPUNIT_ASSERT(_mymaterial);
    _mydata = NULL;

    GenericComponent::name("TestIsotropicLinearGenMaxwellPlaneStrain");

    _mymaterial->PyreComponent::identifier("TestIsotropicLinearGenMaxwellPlaneStrain");
    const char* journal = _mymaterial->PyreComponent::name();
    journal::debug_t debug(journal);
    //debug.activate(); // DEBUGGING
} // setUp


// ----------------------------------------------------------------------
// Deallocate testing data.
void
pylith::materials::TestIsotropicLinearGenMaxwellPlaneStrain::tearDown(void) {
    const char* journal = _mymaterial->PyreComponent::name();
    journal::debug_t debug(journal);
    debug.deactivate(); // DEBUGGING

    TestMaterial::tearDown();

    delete _mymaterial; _mymaterial = NULL;
    delete _mydata; _mydata = NULL;
} // tearDown


// ----------------------------------------------------------------------
// Test useInertia(), useBodyForce(), useReferenceState().
void
pylith::materials::TestIsotropicLinearGenMaxwellPlaneStrain::testAccessors(void) {
    PYLITH_METHOD_BEGIN;

    CPPUNIT_ASSERT(_mymaterial);

    const bool flag = false;

    _mymaterial->useInertia(flag);
	CPPUNIT_ASSERT_EQUAL_MESSAGE("Test of useInertia() failed.", flag, _mymaterial->_useInertia);

    _mymaterial->useInertia(!flag);
    CPPUNIT_ASSERT_EQUAL_MESSAGE("Test of useInertia() failed.", !flag, _mymaterial->_useInertia);

    _mymaterial->useBodyForce(flag);
    CPPUNIT_ASSERT_EQUAL_MESSAGE("Test of useBodyForce() failed.", flag, _mymaterial->_useBodyForce);

    _mymaterial->useBodyForce(!flag);
    CPPUNIT_ASSERT_EQUAL_MESSAGE("Test of useBodyForce() failed.", !flag, _mymaterial->_useBodyForce);

    _mymaterial->useReferenceState(flag);
    CPPUNIT_ASSERT_EQUAL_MESSAGE("Test of useReferenceState() failed.", flag, _mymaterial->_useReferenceState);

    _mymaterial->useReferenceState(!flag);
    CPPUNIT_ASSERT_EQUAL_MESSAGE("Test of useReferenceState() failed.", !flag, _mymaterial->_useReferenceState);

    PYLITH_METHOD_END;
} // testAccessors


// ----------------------------------------------------------------------
// Test auxFieldSetup().
void
pylith::materials::TestIsotropicLinearGenMaxwellPlaneStrain::test_auxFieldSetup(void) {
    PYLITH_METHOD_BEGIN;

    CPPUNIT_ASSERT(_mymaterial);
    CPPUNIT_ASSERT(_mesh);
    CPPUNIT_ASSERT(_mydata);
    CPPUNIT_ASSERT(_mydata->normalizer);
    const PylithReal densityScale = _mydata->normalizer->densityScale();
    const PylithReal lengthScale = _mydata->normalizer->lengthScale();
    const PylithReal timeScale = _mydata->normalizer->timeScale();
    const PylithReal pressureScale = _mydata->normalizer->pressureScale();
    const PylithReal forceScale = pressureScale / lengthScale;
    const PylithReal accelerationScale = lengthScale/(timeScale * timeScale);

    delete _mymaterial->_auxField; _mymaterial->_auxField = new topology::Field(*_mesh); CPPUNIT_ASSERT(_mymaterial->_auxField);
    _mymaterial->_auxFieldSetup();

    // Check discretizations
    { // density
        const char* label = "density";
        const pylith::topology::Field::SubfieldInfo& info = _mymaterial->_auxField->subfieldInfo(label);
        CPPUNIT_ASSERT_EQUAL(size_t(1), info.description.numComponents);
        CPPUNIT_ASSERT_EQUAL(std::string(label), info.description.label);
        CPPUNIT_ASSERT_EQUAL(pylith::topology::Field::SCALAR, info.description.vectorFieldType);
        CPPUNIT_ASSERT_EQUAL(densityScale, info.description.scale);
        CPPUNIT_ASSERT_EQUAL(-1, info.fe.basisOrder);
        CPPUNIT_ASSERT_EQUAL(-1, info.fe.quadOrder);
        CPPUNIT_ASSERT_EQUAL(true, info.fe.isBasisContinuous);
        CPPUNIT_ASSERT_EQUAL(pylith::topology::Field::POLYNOMIAL_SPACE, info.fe.feSpace);
    } // density

    { // shear modulus
        const char* label = "shear_modulus";
        const pylith::topology::Field::SubfieldInfo& info = _mymaterial->_auxField->subfieldInfo(label);
        CPPUNIT_ASSERT_EQUAL(size_t(1), info.description.numComponents);
        CPPUNIT_ASSERT_EQUAL(std::string(label), info.description.label);
        CPPUNIT_ASSERT_EQUAL(pylith::topology::Field::SCALAR, info.description.vectorFieldType);
        CPPUNIT_ASSERT_EQUAL(pressureScale, info.description.scale);
        CPPUNIT_ASSERT_EQUAL(-1, info.fe.basisOrder);
        CPPUNIT_ASSERT_EQUAL(-1, info.fe.quadOrder);
        CPPUNIT_ASSERT_EQUAL(true, info.fe.isBasisContinuous);
        CPPUNIT_ASSERT_EQUAL(pylith::topology::Field::POLYNOMIAL_SPACE, info.fe.feSpace);
    } // shear modulus

    { // bulk modulus
        const char* label = "bulk_modulus";
        const pylith::topology::Field::SubfieldInfo& info = _mymaterial->_auxField->subfieldInfo(label);
        CPPUNIT_ASSERT_EQUAL(size_t(1), info.description.numComponents);
        CPPUNIT_ASSERT_EQUAL(std::string(label), info.description.label);
        CPPUNIT_ASSERT_EQUAL(pylith::topology::Field::SCALAR, info.description.vectorFieldType);
        CPPUNIT_ASSERT_EQUAL(pressureScale, info.description.scale);
        CPPUNIT_ASSERT_EQUAL(-1, info.fe.basisOrder);
        CPPUNIT_ASSERT_EQUAL(-1, info.fe.quadOrder);
        CPPUNIT_ASSERT_EQUAL(true, info.fe.isBasisContinuous);
        CPPUNIT_ASSERT_EQUAL(pylith::topology::Field::POLYNOMIAL_SPACE, info.fe.feSpace);
    } // bulk modulus

    { // Maxwell time 1
        const char* label = "maxwell_time_1";
        const pylith::topology::Field::SubfieldInfo& info = _mymaterial->_auxField->subfieldInfo(label);
        CPPUNIT_ASSERT_EQUAL(size_t(1), info.description.numComponents);
        CPPUNIT_ASSERT_EQUAL(std::string(label), info.description.label);
        CPPUNIT_ASSERT_EQUAL(pylith::topology::Field::SCALAR, info.description.vectorFieldType);
        CPPUNIT_ASSERT_EQUAL(timeScale, info.description.scale);
        CPPUNIT_ASSERT_EQUAL(-1, info.fe.basisOrder);
        CPPUNIT_ASSERT_EQUAL(-1, info.fe.quadOrder);
        CPPUNIT_ASSERT_EQUAL(true, info.fe.isBasisContinuous);
        CPPUNIT_ASSERT_EQUAL(pylith::topology::Field::POLYNOMIAL_SPACE, info.fe.feSpace);
    } // Maxwell time 1

    { // Maxwell time 2
        const char* label = "maxwell_time_2";
        const pylith::topology::Field::SubfieldInfo& info = _mymaterial->_auxField->subfieldInfo(label);
        CPPUNIT_ASSERT_EQUAL(size_t(1), info.description.numComponents);
        CPPUNIT_ASSERT_EQUAL(std::string(label), info.description.label);
        CPPUNIT_ASSERT_EQUAL(pylith::topology::Field::SCALAR, info.description.vectorFieldType);
        CPPUNIT_ASSERT_EQUAL(timeScale, info.description.scale);
        CPPUNIT_ASSERT_EQUAL(-1, info.fe.basisOrder);
        CPPUNIT_ASSERT_EQUAL(-1, info.fe.quadOrder);
        CPPUNIT_ASSERT_EQUAL(true, info.fe.isBasisContinuous);
        CPPUNIT_ASSERT_EQUAL(pylith::topology::Field::POLYNOMIAL_SPACE, info.fe.feSpace);
    } // Maxwell time 2

    { // Maxwell time 3
        const char* label = "maxwell_time_3";
        const pylith::topology::Field::SubfieldInfo& info = _mymaterial->_auxField->subfieldInfo(label);
        CPPUNIT_ASSERT_EQUAL(size_t(1), info.description.numComponents);
        CPPUNIT_ASSERT_EQUAL(std::string(label), info.description.label);
        CPPUNIT_ASSERT_EQUAL(pylith::topology::Field::SCALAR, info.description.vectorFieldType);
        CPPUNIT_ASSERT_EQUAL(timeScale, info.description.scale);
        CPPUNIT_ASSERT_EQUAL(-1, info.fe.basisOrder);
        CPPUNIT_ASSERT_EQUAL(-1, info.fe.quadOrder);
        CPPUNIT_ASSERT_EQUAL(true, info.fe.isBasisContinuous);
        CPPUNIT_ASSERT_EQUAL(pylith::topology::Field::POLYNOMIAL_SPACE, info.fe.feSpace);
    } // Maxwell time 3

    { // Shear modulus ratio 1
        const char* label = "shear_modulus_ratio_1";
        const pylith::topology::Field::SubfieldInfo& info = _mymaterial->_auxField->subfieldInfo(label);
        CPPUNIT_ASSERT_EQUAL(size_t(1), info.description.numComponents);
        CPPUNIT_ASSERT_EQUAL(std::string(label), info.description.label);
        CPPUNIT_ASSERT_EQUAL(pylith::topology::Field::SCALAR, info.description.vectorFieldType);
        CPPUNIT_ASSERT_EQUAL(1.0, info.description.scale);
        CPPUNIT_ASSERT_EQUAL(-1, info.fe.basisOrder);
        CPPUNIT_ASSERT_EQUAL(-1, info.fe.quadOrder);
        CPPUNIT_ASSERT_EQUAL(true, info.fe.isBasisContinuous);
        CPPUNIT_ASSERT_EQUAL(pylith::topology::Field::POLYNOMIAL_SPACE, info.fe.feSpace);
    } // Shear modulus ratio 1

    { // Shear modulus ratio 2
        const char* label = "shear_modulus_ratio_2";
        const pylith::topology::Field::SubfieldInfo& info = _mymaterial->_auxField->subfieldInfo(label);
        CPPUNIT_ASSERT_EQUAL(size_t(1), info.description.numComponents);
        CPPUNIT_ASSERT_EQUAL(std::string(label), info.description.label);
        CPPUNIT_ASSERT_EQUAL(pylith::topology::Field::SCALAR, info.description.vectorFieldType);
        CPPUNIT_ASSERT_EQUAL(1.0, info.description.scale);
        CPPUNIT_ASSERT_EQUAL(-1, info.fe.basisOrder);
        CPPUNIT_ASSERT_EQUAL(-1, info.fe.quadOrder);
        CPPUNIT_ASSERT_EQUAL(true, info.fe.isBasisContinuous);
        CPPUNIT_ASSERT_EQUAL(pylith::topology::Field::POLYNOMIAL_SPACE, info.fe.feSpace);
    } // Shear modulus ratio 2

    { // Shear modulus ratio 3
        const char* label = "shear_modulus_ratio_3";
        const pylith::topology::Field::SubfieldInfo& info = _mymaterial->_auxField->subfieldInfo(label);
        CPPUNIT_ASSERT_EQUAL(size_t(1), info.description.numComponents);
        CPPUNIT_ASSERT_EQUAL(std::string(label), info.description.label);
        CPPUNIT_ASSERT_EQUAL(pylith::topology::Field::SCALAR, info.description.vectorFieldType);
        CPPUNIT_ASSERT_EQUAL(1.0, info.description.scale);
        CPPUNIT_ASSERT_EQUAL(-1, info.fe.basisOrder);
        CPPUNIT_ASSERT_EQUAL(-1, info.fe.quadOrder);
        CPPUNIT_ASSERT_EQUAL(true, info.fe.isBasisContinuous);
        CPPUNIT_ASSERT_EQUAL(pylith::topology::Field::POLYNOMIAL_SPACE, info.fe.feSpace);
    } // Shear modulus ratio 3

    { // Total strain
        const char* label = "total_strain";
        const pylith::topology::Field::SubfieldInfo& info = _mymaterial->_auxField->subfieldInfo(label);
        CPPUNIT_ASSERT_EQUAL(size_t(4), info.description.numComponents);
        CPPUNIT_ASSERT_EQUAL(std::string(label), info.description.label);
        CPPUNIT_ASSERT_EQUAL(pylith::topology::Field::OTHER, info.description.vectorFieldType);
        CPPUNIT_ASSERT_EQUAL(1.0, info.description.scale);
        CPPUNIT_ASSERT_EQUAL(-1, info.fe.basisOrder);
        CPPUNIT_ASSERT_EQUAL(-1, info.fe.quadOrder);
        CPPUNIT_ASSERT_EQUAL(true, info.fe.isBasisContinuous);
        CPPUNIT_ASSERT_EQUAL(pylith::topology::Field::POLYNOMIAL_SPACE, info.fe.feSpace);
    } // Total strain

    { // Viscous strain 1
        const char* label = "viscous_strain_1";
        const pylith::topology::Field::SubfieldInfo& info = _mymaterial->_auxField->subfieldInfo(label);
        CPPUNIT_ASSERT_EQUAL(size_t(4), info.description.numComponents);
        CPPUNIT_ASSERT_EQUAL(std::string(label), info.description.label);
        CPPUNIT_ASSERT_EQUAL(pylith::topology::Field::OTHER, info.description.vectorFieldType);
        CPPUNIT_ASSERT_EQUAL(1.0, info.description.scale);
        CPPUNIT_ASSERT_EQUAL(-1, info.fe.basisOrder);
        CPPUNIT_ASSERT_EQUAL(-1, info.fe.quadOrder);
        CPPUNIT_ASSERT_EQUAL(true, info.fe.isBasisContinuous);
        CPPUNIT_ASSERT_EQUAL(pylith::topology::Field::POLYNOMIAL_SPACE, info.fe.feSpace);
    } // Viscous strain 1

    { // Viscous strain 2
        const char* label = "viscous_strain_2";
        const pylith::topology::Field::SubfieldInfo& info = _mymaterial->_auxField->subfieldInfo(label);
        CPPUNIT_ASSERT_EQUAL(size_t(4), info.description.numComponents);
        CPPUNIT_ASSERT_EQUAL(std::string(label), info.description.label);
        CPPUNIT_ASSERT_EQUAL(pylith::topology::Field::OTHER, info.description.vectorFieldType);
        CPPUNIT_ASSERT_EQUAL(1.0, info.description.scale);
        CPPUNIT_ASSERT_EQUAL(-1, info.fe.basisOrder);
        CPPUNIT_ASSERT_EQUAL(-1, info.fe.quadOrder);
        CPPUNIT_ASSERT_EQUAL(true, info.fe.isBasisContinuous);
        CPPUNIT_ASSERT_EQUAL(pylith::topology::Field::POLYNOMIAL_SPACE, info.fe.feSpace);
    } // Viscous strain 2

    { // Viscous strain 3
        const char* label = "viscous_strain_3";
        const pylith::topology::Field::SubfieldInfo& info = _mymaterial->_auxField->subfieldInfo(label);
        CPPUNIT_ASSERT_EQUAL(size_t(4), info.description.numComponents);
        CPPUNIT_ASSERT_EQUAL(std::string(label), info.description.label);
        CPPUNIT_ASSERT_EQUAL(pylith::topology::Field::OTHER, info.description.vectorFieldType);
        CPPUNIT_ASSERT_EQUAL(1.0, info.description.scale);
        CPPUNIT_ASSERT_EQUAL(-1, info.fe.basisOrder);
        CPPUNIT_ASSERT_EQUAL(-1, info.fe.quadOrder);
        CPPUNIT_ASSERT_EQUAL(true, info.fe.isBasisContinuous);
        CPPUNIT_ASSERT_EQUAL(pylith::topology::Field::POLYNOMIAL_SPACE, info.fe.feSpace);
    } // Viscous strain 3

    if (_mymaterial->_gravityField) { // gravity field
        const char* label = "gravity_field";
        const pylith::topology::Field::SubfieldInfo& info = _mymaterial->_auxField->subfieldInfo(label);
        CPPUNIT_ASSERT_EQUAL(size_t(2), info.description.numComponents);
        CPPUNIT_ASSERT_EQUAL(std::string(label), info.description.label);
        CPPUNIT_ASSERT_EQUAL(pylith::topology::Field::VECTOR, info.description.vectorFieldType);
        CPPUNIT_ASSERT_EQUAL(accelerationScale, info.description.scale);
        CPPUNIT_ASSERT_EQUAL(-1, info.fe.basisOrder);
        CPPUNIT_ASSERT_EQUAL(-1, info.fe.quadOrder);
        CPPUNIT_ASSERT_EQUAL(true, info.fe.isBasisContinuous);
        CPPUNIT_ASSERT_EQUAL(pylith::topology::Field::POLYNOMIAL_SPACE, info.fe.feSpace);
    } // gravity field

    if (_mymaterial->_useBodyForce) { // body force
        const char* label = "body_force";
        const pylith::topology::Field::SubfieldInfo& info = _mymaterial->_auxField->subfieldInfo(label);
        CPPUNIT_ASSERT_EQUAL(size_t(2), info.description.numComponents);
        CPPUNIT_ASSERT_EQUAL(std::string(label), info.description.label);
        CPPUNIT_ASSERT_EQUAL(pylith::topology::Field::VECTOR, info.description.vectorFieldType);
        CPPUNIT_ASSERT_EQUAL(forceScale, info.description.scale);
        CPPUNIT_ASSERT_EQUAL(-1, info.fe.basisOrder);
        CPPUNIT_ASSERT_EQUAL(-1, info.fe.quadOrder);
        CPPUNIT_ASSERT_EQUAL(true, info.fe.isBasisContinuous);
        CPPUNIT_ASSERT_EQUAL(pylith::topology::Field::POLYNOMIAL_SPACE, info.fe.feSpace);
    } // body force

    if (_mymaterial->_useReferenceState) { // reference stress
        const char* label = "reference_stress";
        const pylith::topology::Field::SubfieldInfo& info = _mymaterial->_auxField->subfieldInfo(label);
        CPPUNIT_ASSERT_EQUAL(size_t(4), info.description.numComponents);
        CPPUNIT_ASSERT_EQUAL(std::string(label), info.description.label);
        CPPUNIT_ASSERT_EQUAL(pylith::topology::Field::OTHER, info.description.vectorFieldType);
        CPPUNIT_ASSERT_EQUAL(pressureScale, info.description.scale);
        CPPUNIT_ASSERT_EQUAL(-1, info.fe.basisOrder);
        CPPUNIT_ASSERT_EQUAL(-1, info.fe.quadOrder);
        CPPUNIT_ASSERT_EQUAL(true, info.fe.isBasisContinuous);
        CPPUNIT_ASSERT_EQUAL(pylith::topology::Field::POLYNOMIAL_SPACE, info.fe.feSpace);
    } // reference stress

    if (_mymaterial->_useReferenceState) { // reference strain
        const char* label = "reference_strain";
        const pylith::topology::Field::SubfieldInfo& info = _mymaterial->_auxField->subfieldInfo(label);
        CPPUNIT_ASSERT_EQUAL(size_t(4), info.description.numComponents);
        CPPUNIT_ASSERT_EQUAL(std::string(label), info.description.label);
        CPPUNIT_ASSERT_EQUAL(pylith::topology::Field::OTHER, info.description.vectorFieldType);
        CPPUNIT_ASSERT_EQUAL(1.0, info.description.scale);
        CPPUNIT_ASSERT_EQUAL(-1, info.fe.basisOrder);
        CPPUNIT_ASSERT_EQUAL(-1, info.fe.quadOrder);
        CPPUNIT_ASSERT_EQUAL(true, info.fe.isBasisContinuous);
        CPPUNIT_ASSERT_EQUAL(pylith::topology::Field::POLYNOMIAL_SPACE, info.fe.feSpace);
    } // reference strain

    PYLITH_METHOD_END;
} // test_auxFieldSetup


// ----------------------------------------------------------------------
// Test getAuxField().
void
pylith::materials::TestIsotropicLinearGenMaxwellPlaneStrain::testGetAuxField(void) {
    PYLITH_METHOD_BEGIN;

    _initializeFull();

    CPPUNIT_ASSERT(_mymaterial);
    CPPUNIT_ASSERT(_mesh);
    CPPUNIT_ASSERT(_mydata->normalizer);
    const PylithReal lengthScale = _mydata->normalizer->lengthScale();

    const pylith::topology::Field& auxField = _mymaterial->auxField();
    { // Test getting density field.
        pylith::topology::Field density(*_mesh);
        density.copySubfield(auxField, "density");

        //density.view("DENSITY"); // DEBUGGING

        // Check result
        CPPUNIT_ASSERT_EQUAL(std::string("density"), std::string(density.label()));
        CPPUNIT_ASSERT_EQUAL(_mydata->dimension, density.spaceDim());

        pylith::topology::FieldQuery queryDensity(density);
        queryDensity.initializeWithDefaultQueryFns();
        queryDensity.openDB(_mydata->auxDB, lengthScale);

        PylithReal norm = 0.0;
        const PylithReal t = _mydata->t;
        const PetscDM dm = density.dmMesh(); CPPUNIT_ASSERT(dm);
        PetscErrorCode err = DMPlexComputeL2DiffLocal(dm, t, queryDensity.functions(), (void**)queryDensity.contextPtrs(), density.localVector(), &norm); CPPUNIT_ASSERT(!err);
        queryDensity.closeDB(_mydata->auxDB);

        const PylithReal tolerance = 1.0e-6;
		CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("Test extracting density subfield from auxiliary field failed.", 0.0, norm, tolerance);
    } // Test getting density field

    { // Test getting bulk_modulus field.
        pylith::topology::Field bulkModulus(*_mesh);
        bulkModulus.copySubfield(auxField, "bulk_modulus");

        //bulkModulus.view("BULK MODULUS"); // DEBUGGING

        // Check result
        CPPUNIT_ASSERT_EQUAL(std::string("bulk_modulus"), std::string(bulkModulus.label()));
        CPPUNIT_ASSERT_EQUAL(_mydata->dimension, bulkModulus.spaceDim());

        pylith::topology::FieldQuery queryBulkModulus(bulkModulus);
        queryBulkModulus.initializeWithDefaultQueryFns();
        queryBulkModulus.openDB(_mydata->auxDB, lengthScale);

        PylithReal norm = 0.0;
        const PylithReal t = _mydata->t;
        const PetscDM dm = bulkModulus.dmMesh(); CPPUNIT_ASSERT(dm);
        PetscErrorCode err = DMPlexComputeL2DiffLocal(dm, t, queryBulkModulus.functions(), (void**)queryBulkModulus.contextPtrs(), bulkModulus.localVector(), &norm); CPPUNIT_ASSERT(!err);
        queryBulkModulus.closeDB(_mydata->auxDB);

        const PylithReal tolerance = 1.0e-6;
		CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("Test extracting bulk modulus subfield from auxiliary field failed.", 0.0, norm, tolerance);
    } // Test getting bulkModulus field

    { // Test getting maxwellTime1 field.
        pylith::topology::Field maxwellTime1(*_mesh);
        maxwellTime1.copySubfield(auxField, "maxwell_time_1");

        //maxwellTime1.view("MAXWELL TIME 1"); // DEBUGGING

        // Check result
        CPPUNIT_ASSERT_EQUAL(std::string("maxwell_time_1"), std::string(maxwellTime1.label()));
        CPPUNIT_ASSERT_EQUAL(_mydata->dimension, maxwellTime1.spaceDim());

        pylith::topology::FieldQuery queryMaxwellTime1(maxwellTime1);
        queryMaxwellTime1.initializeWithDefaultQueryFns();
        queryMaxwellTime1.openDB(_mydata->auxDB, lengthScale);

        PylithReal norm = 0.0;
        const PylithReal t = _mydata->t;
        const PetscDM dm = maxwellTime1.dmMesh(); CPPUNIT_ASSERT(dm);
        PetscErrorCode err = DMPlexComputeL2DiffLocal(dm, t, queryMaxwellTime1.functions(), (void**)queryMaxwellTime1.contextPtrs(), maxwellTime1.localVector(), &norm); CPPUNIT_ASSERT(!err);
        queryMaxwellTime1.closeDB(_mydata->auxDB);

        const PylithReal tolerance = 1.0e-6;
		CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("Test extracting Maxwell time 1 subfield from auxiliary field failed.", 0.0, norm, tolerance);
    } // Test getting maxwellTime1 field

    { // Test getting maxwellTime2 field.
        pylith::topology::Field maxwellTime2(*_mesh);
        maxwellTime2.copySubfield(auxField, "maxwell_time_2");

        //maxwellTime2.view("MAXWELL TIME 2"); // DEBUGGING

        // Check result
        CPPUNIT_ASSERT_EQUAL(std::string("maxwell_time_2"), std::string(maxwellTime2.label()));
        CPPUNIT_ASSERT_EQUAL(_mydata->dimension, maxwellTime2.spaceDim());

        pylith::topology::FieldQuery queryMaxwellTime2(maxwellTime2);
        queryMaxwellTime2.initializeWithDefaultQueryFns();
        queryMaxwellTime2.openDB(_mydata->auxDB, lengthScale);

        PylithReal norm = 0.0;
        const PylithReal t = _mydata->t;
        const PetscDM dm = maxwellTime2.dmMesh(); CPPUNIT_ASSERT(dm);
        PetscErrorCode err = DMPlexComputeL2DiffLocal(dm, t, queryMaxwellTime2.functions(), (void**)queryMaxwellTime2.contextPtrs(), maxwellTime2.localVector(), &norm); CPPUNIT_ASSERT(!err);
        queryMaxwellTime2.closeDB(_mydata->auxDB);

        const PylithReal tolerance = 1.0e-6;
		CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("Test extracting Maxwell time 2 subfield from auxiliary field failed.", 0.0, norm, tolerance);
    } // Test getting maxwellTime2 field

    { // Test getting maxwellTime3 field.
        pylith::topology::Field maxwellTime3(*_mesh);
        maxwellTime3.copySubfield(auxField, "maxwell_time_3");

        //maxwellTime3.view("MAXWELL TIME 3"); // DEBUGGING

        // Check result
        CPPUNIT_ASSERT_EQUAL(std::string("maxwell_time_3"), std::string(maxwellTime3.label()));
        CPPUNIT_ASSERT_EQUAL(_mydata->dimension, maxwellTime3.spaceDim());

        pylith::topology::FieldQuery queryMaxwellTime3(maxwellTime3);
        queryMaxwellTime3.initializeWithDefaultQueryFns();
        queryMaxwellTime3.openDB(_mydata->auxDB, lengthScale);

        PylithReal norm = 0.0;
        const PylithReal t = _mydata->t;
        const PetscDM dm = maxwellTime3.dmMesh(); CPPUNIT_ASSERT(dm);
        PetscErrorCode err = DMPlexComputeL2DiffLocal(dm, t, queryMaxwellTime3.functions(), (void**)queryMaxwellTime3.contextPtrs(), maxwellTime3.localVector(), &norm); CPPUNIT_ASSERT(!err);
        queryMaxwellTime3.closeDB(_mydata->auxDB);

        const PylithReal tolerance = 1.0e-6;
		CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("Test extracting Maxwell time 3 subfield from auxiliary field failed.", 0.0, norm, tolerance);
    } // Test getting maxwellTime3 field

    { // Test getting shearModulusRatio1 field.
        pylith::topology::Field shearModulusRatio1(*_mesh);
        shearModulusRatio1.copySubfield(auxField, "shear_modulus_ratio_1");

        //shearModulusRatio1.view("SHEAR MODULUS RATIO 1"); // DEBUGGING

        // Check result
        CPPUNIT_ASSERT_EQUAL(std::string("shear_modulus_ratio_1"), std::string(shearModulusRatio1.label()));
        CPPUNIT_ASSERT_EQUAL(_mydata->dimension, shearModulusRatio1.spaceDim());

        pylith::topology::FieldQuery queryShearModulusRatio1(shearModulusRatio1);
        queryShearModulusRatio1.initializeWithDefaultQueryFns();
        queryShearModulusRatio1.openDB(_mydata->auxDB, lengthScale);

        PylithReal norm = 0.0;
        const PylithReal t = _mydata->t;
        const PetscDM dm = shearModulusRatio1.dmMesh(); CPPUNIT_ASSERT(dm);
        PetscErrorCode err = DMPlexComputeL2DiffLocal(dm, t, queryShearModulusRatio1.functions(), (void**)queryShearModulusRatio1.contextPtrs(), shearModulusRatio1.localVector(), &norm); CPPUNIT_ASSERT(!err);
        queryShearModulusRatio1.closeDB(_mydata->auxDB);

        const PylithReal tolerance = 1.0e-6;
		CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("Test extracting shear modulus ratio 1 subfield from auxiliary field failed.", 0.0, norm, tolerance);
    } // Test getting shearModulusRatio1 field

    { // Test getting shearModulusRatio2 field.
        pylith::topology::Field shearModulusRatio2(*_mesh);
        shearModulusRatio2.copySubfield(auxField, "shear_modulus_ratio_2");

        //shearModulusRatio2.view("SHEAR MODULUS RATIO 2"); // DEBUGGING

        // Check result
        CPPUNIT_ASSERT_EQUAL(std::string("shear_modulus_ratio_2"), std::string(shearModulusRatio2.label()));
        CPPUNIT_ASSERT_EQUAL(_mydata->dimension, shearModulusRatio2.spaceDim());

        pylith::topology::FieldQuery queryShearModulusRatio2(shearModulusRatio2);
        queryShearModulusRatio2.initializeWithDefaultQueryFns();
        queryShearModulusRatio2.openDB(_mydata->auxDB, lengthScale);

        PylithReal norm = 0.0;
        const PylithReal t = _mydata->t;
        const PetscDM dm = shearModulusRatio2.dmMesh(); CPPUNIT_ASSERT(dm);
        PetscErrorCode err = DMPlexComputeL2DiffLocal(dm, t, queryShearModulusRatio2.functions(), (void**)queryShearModulusRatio2.contextPtrs(), shearModulusRatio2.localVector(), &norm); CPPUNIT_ASSERT(!err);
        queryShearModulusRatio2.closeDB(_mydata->auxDB);

        const PylithReal tolerance = 1.0e-6;
		CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("Test extracting shear modulus ratio 2 subfield from auxiliary field failed.", 0.0, norm, tolerance);
    } // Test getting shearModulusRatio2 field

    { // Test getting shearModulusRatio3 field.
        pylith::topology::Field shearModulusRatio3(*_mesh);
        shearModulusRatio3.copySubfield(auxField, "shear_modulus_ratio_3");

        //shearModulusRatio3.view("SHEAR MODULUS RATIO 3"); // DEBUGGING

        // Check result
        CPPUNIT_ASSERT_EQUAL(std::string("shear_modulus_ratio_3"), std::string(shearModulusRatio3.label()));
        CPPUNIT_ASSERT_EQUAL(_mydata->dimension, shearModulusRatio3.spaceDim());

        pylith::topology::FieldQuery queryShearModulusRatio3(shearModulusRatio3);
        queryShearModulusRatio3.initializeWithDefaultQueryFns();
        queryShearModulusRatio3.openDB(_mydata->auxDB, lengthScale);

        PylithReal norm = 0.0;
        const PylithReal t = _mydata->t;
        const PetscDM dm = shearModulusRatio3.dmMesh(); CPPUNIT_ASSERT(dm);
        PetscErrorCode err = DMPlexComputeL2DiffLocal(dm, t, queryShearModulusRatio3.functions(), (void**)queryShearModulusRatio3.contextPtrs(), shearModulusRatio3.localVector(), &norm); CPPUNIT_ASSERT(!err);
        queryShearModulusRatio3.closeDB(_mydata->auxDB);

        const PylithReal tolerance = 1.0e-6;
		CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("Test extracting shear modulus ratio 3 subfield from auxiliary field failed.", 0.0, norm, tolerance);
    } // Test getting shearModulusRatio3 field

	if (_mymaterial->_useReferenceState) { // Test getting reference_stress field.
        pylith::topology::Field referenceStress(*_mesh);
        referenceStress.copySubfield(auxField, "reference_stress");

        //referenceStress.view("REFERENCE STRESS"); // DEBUGGING

        // Check result
        CPPUNIT_ASSERT_EQUAL(std::string("reference_stress"), std::string(referenceStress.label()));
        CPPUNIT_ASSERT_EQUAL(_mydata->dimension, referenceStress.spaceDim());

        pylith::topology::FieldQuery queryRefStress(referenceStress);
        queryRefStress.initializeWithDefaultQueryFns();
        queryRefStress.openDB(_mydata->auxDB, lengthScale);

        PylithReal norm = 0.0;
        const PylithReal t = _mydata->t;
        const PetscDM dm = referenceStress.dmMesh(); CPPUNIT_ASSERT(dm);
        PetscErrorCode err = DMPlexComputeL2DiffLocal(dm, t, queryRefStress.functions(), (void**)queryRefStress.contextPtrs(), referenceStress.localVector(), &norm); CPPUNIT_ASSERT(!err);
        queryRefStress.closeDB(_mydata->auxDB);

        const PylithReal tolerance = 1.0e-6;
        CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("Test extracting reference stress subfield from auxiliary field failed.", 0.0, norm, tolerance);
    } // Test getting reference_stress field

    // PYLITH_JOURNAL_WARNING(":TODO: @charles Add test for getting reference stress/strain.");

    PYLITH_METHOD_END;
} // testGetAuxField


// ----------------------------------------------------------------------
// Get material.
pylith::materials::Material*
pylith::materials::TestIsotropicLinearGenMaxwellPlaneStrain::_material(void)
{ // _material
    return _mymaterial;
} // _material


// ----------------------------------------------------------------------
// Get test data.
pylith::materials::TestMaterial_Data*
pylith::materials::TestIsotropicLinearGenMaxwellPlaneStrain::_data(void)
{ // _data
    return _mydata;
} // _data


// ----------------------------------------------------------------------
// Setup and populate solution fields.
void
pylith::materials::TestIsotropicLinearGenMaxwellPlaneStrain::_setupSolutionFields(void)
{ // _setupSolutionFields
    PYLITH_METHOD_BEGIN;

    CPPUNIT_ASSERT(_solutionFields);

    CPPUNIT_ASSERT( (!_mydata->isExplicit && 1 == _mydata->numSolnSubfields) ||
                    (_mydata->isExplicit && 2 == _mydata->numSolnSubfields) );
    CPPUNIT_ASSERT(_mydata->solnDiscretizations);
    CPPUNIT_ASSERT(_mydata->normalizer);

    { // Solution
        pylith::topology::Field& solution = _solutionFields->get("solution");
        pylith::problems::SolutionFactory factory(solution, *_mydata->normalizer);
        factory.displacement(_mydata->solnDiscretizations[0]);
        if (_mydata->isExplicit) {
            factory.velocity(_mydata->solnDiscretizations[1]);
        } // if
        solution.subfieldsSetup();
        solution.allocate();
        factory.setValues(_mydata->solnDB);
    } // Solution

    { // Time derivative of solution
        pylith::topology::Field& solutionDot = _solutionFields->get("solution_dot");
        pylith::problems::SolutionFactory factory(solutionDot, *_mydata->normalizer);
        factory.displacementDot(_mydata->solnDiscretizations[0]);
        if (_mydata->isExplicit) {
            factory.velocityDot(_mydata->solnDiscretizations[1]);
        } // if
        solutionDot.subfieldsSetup();
        solutionDot.allocate();
        factory.setValues(_mydata->solnDB);
    } // Time derivative of solution

    { // Perturbation
        pylith::topology::Field& perturbation = _solutionFields->get("perturbation");
        const pylith::topology::Field& solution = _solutionFields->get("solution");
        perturbation.cloneSection(solution);
        perturbation.allocate();
        perturbation.zeroLocal();
		pylith::problems::SolutionFactory factory(perturbation, *_mydata->normalizer);
        factory.setValues(_mydata->perturbDB);
    } // Perturbation @ t2

    { // Time derivative of solution @ t2
        pylith::topology::Field& perturbationDot = _solutionFields->get("perturbation_dot");
        const pylith::topology::Field& solutionDot = _solutionFields->get("solution_dot");
        perturbationDot.cloneSection(solutionDot);
        perturbationDot.allocate();
        perturbationDot.zeroLocal();
		pylith::problems::SolutionFactory factory(perturbationDot, *_mydata->normalizer);
        factory.setValues(_mydata->perturbDB);
    } // Time derivative of solution @ t2

    PYLITH_METHOD_END;
} // _setupSolutionFields

// ----------------------------------------------------------------------
// Constructor
pylith::materials::TestIsotropicLinearGenMaxwellPlaneStrain_Data::TestIsotropicLinearGenMaxwellPlaneStrain_Data(void)
{ // constructor
    dimension = 2;
    gravityVector[0] = 0.0; // Use scales in test to provide correct nondimensional value.
    gravityVector[1] = 0.0;
    gravityVector[2] = 0;

    cs = new spatialdata::geocoords::CSCart; CPPUNIT_ASSERT(cs);
    cs->setSpaceDim(dimension);
    cs->initialize();

	// Some auxiliary subfields get updated in updateStateVars().
    auxUpdateDB = new spatialdata::spatialdb::UserFunctionDB; CPPUNIT_ASSERT(auxUpdateDB);

    solnDB->coordsys(*cs);
    perturbDB->coordsys(*cs);
    auxDB->coordsys(*cs);
    auxUpdateDB->coordsys(*cs);
} // constructor


// ----------------------------------------------------------------------
// Destructor
pylith::materials::TestIsotropicLinearGenMaxwellPlaneStrain_Data::~TestIsotropicLinearGenMaxwellPlaneStrain_Data(void) {}


// End of file
