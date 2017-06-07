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

#include "TestIsotropicLinearMaxwellPlaneStrain.hh" // Implementation of class methods

#include "pylith/materials/IsotropicLinearMaxwellPlaneStrain.hh" // USES IsotropicLinearMaxwellPlaneStrain
#include "pylith/materials/Query.hh" // USES Query

#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/Fields.hh" // USES Fields
#include "pylith/topology/FieldQuery.hh" // USES FieldQuery
#include "pylith/problems/SolutionFactory.hh" // USES SolutionFactory

#include "pylith/utils/error.hh" // USES PYLITH_METHOD_BEGIN/END
#include "pylith/utils/journals.hh" // USES JournalingComponent

#include "spatialdata/spatialdb/SimpleGridDB.hh" // USES SimpleGridDB
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

extern "C" {
#include "pylith/fekernels/dispvel.h" // USES DispVel kernels
#include "pylith/fekernels/elasticity.h" // USES Elasticity kernels
#include "pylith/fekernels/linearmaxwellplanestrain.h" // USES IsotropicLinearMaxwellPlaneStrain kernels
}


// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::materials::TestIsotropicLinearMaxwellPlaneStrain::setUp(void)
{ // setUp
    TestMaterialNew::setUp();
    _mymaterial = new IsotropicLinearMaxwellPlaneStrain(); CPPUNIT_ASSERT(_mymaterial);
    _mydata = NULL;

    _mymaterial->PyreComponent::identifier("TestIsotropicLinearMaxwellPlaneStrain");
    const char* journal = _mymaterial->PyreComponent::name();
    journal::debug_t debug(journal);
    //debug.activate(); // DEBUGGING
} // setUp


// ----------------------------------------------------------------------
// Deallocate testing data.
void
pylith::materials::TestIsotropicLinearMaxwellPlaneStrain::tearDown(void)
{ // tearDown
    const char* journal = _mymaterial->PyreComponent::name();
    journal::debug_t debug(journal);
    debug.deactivate(); // DEBUGGING

    TestMaterialNew::tearDown();

    delete _mymaterial; _mymaterial = NULL;
    delete _mydata; _mydata = NULL;
} // tearDown


// ----------------------------------------------------------------------
// Test useInertia().
void
pylith::materials::TestIsotropicLinearMaxwellPlaneStrain::testUseInertia(void)
{ // testUseInertia
    PYLITH_METHOD_BEGIN;

    CPPUNIT_ASSERT(_mymaterial);

    const bool flag = false; // default
    CPPUNIT_ASSERT_EQUAL(flag, _mymaterial->_useInertia);

    _mymaterial->useInertia(!flag);
    CPPUNIT_ASSERT_EQUAL(!flag, _mymaterial->_useInertia);

    PYLITH_METHOD_END;
} // testUseInertia


// ----------------------------------------------------------------------
// Test useBodyForce().
void
pylith::materials::TestIsotropicLinearMaxwellPlaneStrain::testUseBodyForce(void)
{ // testUseBodyForce
    PYLITH_METHOD_BEGIN;

    CPPUNIT_ASSERT(_mymaterial);

    const bool flag = false; // default
    CPPUNIT_ASSERT_EQUAL(flag, _mymaterial->_useBodyForce);

    _mymaterial->useBodyForce(!flag);
    CPPUNIT_ASSERT_EQUAL(!flag, _mymaterial->_useBodyForce);

    PYLITH_METHOD_END;
} // testUseBodyForce


// ----------------------------------------------------------------------
// Test useReferenceState().
void
pylith::materials::TestIsotropicLinearMaxwellPlaneStrain::testUseReferenceState(void)
{ // testUseReferenceState
    PYLITH_METHOD_BEGIN;

    CPPUNIT_ASSERT(_mymaterial);

    const bool flag = false; // default
    CPPUNIT_ASSERT_EQUAL(flag, _mymaterial->_useReferenceState);

    _mymaterial->useReferenceState(!flag);
    CPPUNIT_ASSERT_EQUAL(!flag, _mymaterial->_useReferenceState);

    PYLITH_METHOD_END;
} // testUseReferenceState


// ----------------------------------------------------------------------
// Test auxFieldsSetup().
void
pylith::materials::TestIsotropicLinearMaxwellPlaneStrain::test_auxFieldsSetup(void)
{ // test_auxFieldsSetup
    PYLITH_METHOD_BEGIN;

    CPPUNIT_ASSERT(_mymaterial);
    CPPUNIT_ASSERT(_mesh);
    CPPUNIT_ASSERT(_mydata);
    CPPUNIT_ASSERT(_mydata->normalizer);
    const PylithReal densityScale = _mydata->normalizer->densityScale();
    const PylithReal lengthScale = _mydata->normalizer->lengthScale();
    const PylithReal timeScale = _mydata->normalizer->timeScale();
    const PylithReal pressureScale = _mydata->normalizer->pressureScale();
    const PylithReal forceScale = densityScale * lengthScale / (timeScale * timeScale);
    const PylithReal accelerationScale = lengthScale/(timeScale * timeScale);

    delete _mymaterial->_auxFields; _mymaterial->_auxFields = new topology::Field(*_mesh); CPPUNIT_ASSERT(_mymaterial->_auxFields);
    delete _mymaterial->_auxFieldsQuery; _mymaterial->_auxFieldsQuery = new topology::FieldQuery(*_mymaterial->_auxFields); CPPUNIT_ASSERT(_mymaterial->_auxFieldsQuery);
    _mymaterial->_auxFieldsSetup();

    // Check discretizations
    { // density
        const char* label = "density";
        const pylith::topology::Field::SubfieldInfo& info = _mymaterial->_auxFields->subfieldInfo(label);
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
        const pylith::topology::Field::SubfieldInfo& info = _mymaterial->_auxFields->subfieldInfo(label);
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
        const pylith::topology::Field::SubfieldInfo& info = _mymaterial->_auxFields->subfieldInfo(label);
        CPPUNIT_ASSERT_EQUAL(size_t(1), info.description.numComponents);
        CPPUNIT_ASSERT_EQUAL(std::string(label), info.description.label);
        CPPUNIT_ASSERT_EQUAL(pylith::topology::Field::SCALAR, info.description.vectorFieldType);
        CPPUNIT_ASSERT_EQUAL(pressureScale, info.description.scale);
        CPPUNIT_ASSERT_EQUAL(-1, info.fe.basisOrder);
        CPPUNIT_ASSERT_EQUAL(-1, info.fe.quadOrder);
        CPPUNIT_ASSERT_EQUAL(true, info.fe.isBasisContinuous);
        CPPUNIT_ASSERT_EQUAL(pylith::topology::Field::POLYNOMIAL_SPACE, info.fe.feSpace);
    } // bulk modulus

    { // Maxwell time
        const char* label = "maxwell_time";
        const pylith::topology::Field::SubfieldInfo& info = _mymaterial->_auxFields->subfieldInfo(label);
        CPPUNIT_ASSERT_EQUAL(size_t(1), info.description.numComponents);
        CPPUNIT_ASSERT_EQUAL(std::string(label), info.description.label);
        CPPUNIT_ASSERT_EQUAL(pylith::topology::Field::SCALAR, info.description.vectorFieldType);
        CPPUNIT_ASSERT_EQUAL(timeScale, info.description.scale);
        CPPUNIT_ASSERT_EQUAL(-1, info.fe.basisOrder);
        CPPUNIT_ASSERT_EQUAL(-1, info.fe.quadOrder);
        CPPUNIT_ASSERT_EQUAL(true, info.fe.isBasisContinuous);
        CPPUNIT_ASSERT_EQUAL(pylith::topology::Field::POLYNOMIAL_SPACE, info.fe.feSpace);
    } // Maxwell time

    { // Total strain
        const char* label = "total_strain";
        const pylith::topology::Field::SubfieldInfo& info = _mymaterial->_auxFields->subfieldInfo(label);
        CPPUNIT_ASSERT_EQUAL(size_t(4), info.description.numComponents);
        CPPUNIT_ASSERT_EQUAL(std::string(label), info.description.label);
        CPPUNIT_ASSERT_EQUAL(pylith::topology::Field::OTHER, info.description.vectorFieldType);
        CPPUNIT_ASSERT_EQUAL(1.0, info.description.scale);
        CPPUNIT_ASSERT_EQUAL(-1, info.fe.basisOrder);
        CPPUNIT_ASSERT_EQUAL(-1, info.fe.quadOrder);
        CPPUNIT_ASSERT_EQUAL(true, info.fe.isBasisContinuous);
        CPPUNIT_ASSERT_EQUAL(pylith::topology::Field::POLYNOMIAL_SPACE, info.fe.feSpace);
    } // Total strain

    { // Viscous strain
        const char* label = "viscous_strain";
        const pylith::topology::Field::SubfieldInfo& info = _mymaterial->_auxFields->subfieldInfo(label);
        CPPUNIT_ASSERT_EQUAL(size_t(4), info.description.numComponents);
        CPPUNIT_ASSERT_EQUAL(std::string(label), info.description.label);
        CPPUNIT_ASSERT_EQUAL(pylith::topology::Field::OTHER, info.description.vectorFieldType);
        CPPUNIT_ASSERT_EQUAL(1.0, info.description.scale);
        CPPUNIT_ASSERT_EQUAL(-1, info.fe.basisOrder);
        CPPUNIT_ASSERT_EQUAL(-1, info.fe.quadOrder);
        CPPUNIT_ASSERT_EQUAL(true, info.fe.isBasisContinuous);
        CPPUNIT_ASSERT_EQUAL(pylith::topology::Field::POLYNOMIAL_SPACE, info.fe.feSpace);
    } // Viscous strain

    if (_mydata->useBodyForce) { // body force
        const char* label = "body_force";
        const pylith::topology::Field::SubfieldInfo& info = _mymaterial->_auxFields->subfieldInfo(label);
        CPPUNIT_ASSERT_EQUAL(size_t(2), info.description.numComponents);
        CPPUNIT_ASSERT_EQUAL(std::string(label), info.description.label);
        CPPUNIT_ASSERT_EQUAL(pylith::topology::Field::VECTOR, info.description.vectorFieldType);
        CPPUNIT_ASSERT_EQUAL(forceScale, info.description.scale);
        CPPUNIT_ASSERT_EQUAL(-1, info.fe.basisOrder);
        CPPUNIT_ASSERT_EQUAL(-1, info.fe.quadOrder);
        CPPUNIT_ASSERT_EQUAL(true, info.fe.isBasisContinuous);
        CPPUNIT_ASSERT_EQUAL(pylith::topology::Field::POLYNOMIAL_SPACE, info.fe.feSpace);
    } // body force

    if (_mydata->useGravity) { // gravity field
        const char* label = "gravity_field";
        const pylith::topology::Field::SubfieldInfo& info = _mymaterial->_auxFields->subfieldInfo(label);
        CPPUNIT_ASSERT_EQUAL(size_t(2), info.description.numComponents);
        CPPUNIT_ASSERT_EQUAL(std::string(label), info.description.label);
        CPPUNIT_ASSERT_EQUAL(pylith::topology::Field::VECTOR, info.description.vectorFieldType);
        CPPUNIT_ASSERT_EQUAL(accelerationScale, info.description.scale);
        CPPUNIT_ASSERT_EQUAL(-1, info.fe.basisOrder);
        CPPUNIT_ASSERT_EQUAL(-1, info.fe.quadOrder);
        CPPUNIT_ASSERT_EQUAL(true, info.fe.isBasisContinuous);
        CPPUNIT_ASSERT_EQUAL(pylith::topology::Field::POLYNOMIAL_SPACE, info.fe.feSpace);
    } // gravity field

    if (_mydata->useReferenceState) { // reference stress
        const char* label = "reference_stress";
        const pylith::topology::Field::SubfieldInfo& info = _mymaterial->_auxFields->subfieldInfo(label);
        CPPUNIT_ASSERT_EQUAL(size_t(4), info.description.numComponents);
        CPPUNIT_ASSERT_EQUAL(std::string(label), info.description.label);
        CPPUNIT_ASSERT_EQUAL(pylith::topology::Field::OTHER, info.description.vectorFieldType);
        CPPUNIT_ASSERT_EQUAL(pressureScale, info.description.scale);
        CPPUNIT_ASSERT_EQUAL(-1, info.fe.basisOrder);
        CPPUNIT_ASSERT_EQUAL(-1, info.fe.quadOrder);
        CPPUNIT_ASSERT_EQUAL(true, info.fe.isBasisContinuous);
        CPPUNIT_ASSERT_EQUAL(pylith::topology::Field::POLYNOMIAL_SPACE, info.fe.feSpace);
    } // reference stress

    if (_mydata->useReferenceState) { // reference strain
        const char* label = "reference_strain";
        const pylith::topology::Field::SubfieldInfo& info = _mymaterial->_auxFields->subfieldInfo(label);
        CPPUNIT_ASSERT_EQUAL(size_t(4), info.description.numComponents);
        CPPUNIT_ASSERT_EQUAL(std::string(label), info.description.label);
        CPPUNIT_ASSERT_EQUAL(pylith::topology::Field::OTHER, info.description.vectorFieldType);
        CPPUNIT_ASSERT_EQUAL(1.0, info.description.scale);
        CPPUNIT_ASSERT_EQUAL(-1, info.fe.basisOrder);
        CPPUNIT_ASSERT_EQUAL(-1, info.fe.quadOrder);
        CPPUNIT_ASSERT_EQUAL(true, info.fe.isBasisContinuous);
        CPPUNIT_ASSERT_EQUAL(pylith::topology::Field::POLYNOMIAL_SPACE, info.fe.feSpace);
    } // reference strain

    // Make sure DB query functions are set correctly.
    CPPUNIT_ASSERT_EQUAL(&pylith::topology::FieldQuery::dbQueryGeneric, _mymaterial->_auxFieldsQuery->queryFn("density"));
    CPPUNIT_ASSERT_EQUAL(&pylith::materials::Query::dbQueryShearModulus, _mymaterial->_auxFieldsQuery->queryFn("shear_modulus"));
    CPPUNIT_ASSERT_EQUAL(&pylith::materials::Query::dbQueryBulkModulus, _mymaterial->_auxFieldsQuery->queryFn("bulk_modulus"));
    CPPUNIT_ASSERT_EQUAL(&pylith::materials::Query::dbQueryMaxwellTime, _mymaterial->_auxFieldsQuery->queryFn("maxwell_time"));
    CPPUNIT_ASSERT_EQUAL(&pylith::topology::FieldQuery::dbQueryGeneric, _mymaterial->_auxFieldsQuery->queryFn("total_strain"));
    CPPUNIT_ASSERT_EQUAL(&pylith::topology::FieldQuery::dbQueryGeneric, _mymaterial->_auxFieldsQuery->queryFn("viscous_strain"));
    if (_mydata->useGravity) {
        CPPUNIT_ASSERT_EQUAL(&pylith::materials::Query::dbQueryGravityField, _mymaterial->_auxFieldsQuery->queryFn("gravity_field"));
    } // if
    if (_mydata->useBodyForce) {
        CPPUNIT_ASSERT_EQUAL(&pylith::topology::FieldQuery::dbQueryGeneric, _mymaterial->_auxFieldsQuery->queryFn("body_force"));
    } // if
    if (_mydata->useReferenceState) {
        CPPUNIT_ASSERT_EQUAL(&pylith::topology::FieldQuery::dbQueryGeneric, _mymaterial->_auxFieldsQuery->queryFn("reference_stress"));
        CPPUNIT_ASSERT_EQUAL(&pylith::topology::FieldQuery::dbQueryGeneric, _mymaterial->_auxFieldsQuery->queryFn("reference_strain"));
    } // if

    PYLITH_METHOD_END;
} // test_auxFieldsSetup


// ----------------------------------------------------------------------
// Test getAuxField().
void
pylith::materials::TestIsotropicLinearMaxwellPlaneStrain::testGetAuxField(void)
{ // testGetAuxField
    PYLITH_METHOD_BEGIN;

    _initializeFull();

    CPPUNIT_ASSERT(_mymaterial);
    CPPUNIT_ASSERT(_mesh);
    CPPUNIT_ASSERT(_mydata->normalizer);
    const PylithReal lengthScale = _mydata->normalizer->lengthScale();

    { // Test getting density field.
        pylith::topology::Field density(*_mesh);
        _mymaterial->getAuxField(&density, "density");

        //density.view("DENSITY"); // DEBUGGING

        // Check result
        CPPUNIT_ASSERT_EQUAL(std::string("density"), std::string(density.label()));
        CPPUNIT_ASSERT_EQUAL(_mydata->dimension, density.spaceDim());

        pylith::topology::FieldQuery queryDensity(density);
        queryDensity.queryFn("density", pylith::topology::FieldQuery::dbQueryGeneric);
        queryDensity.openDB(_auxDB, lengthScale);

        PylithReal norm = 0.0;
        const PylithReal t = _mydata->t;
        const PetscDM dm = density.dmMesh(); CPPUNIT_ASSERT(dm);
        PetscErrorCode err = DMComputeL2Diff(dm, t, queryDensity.functions(), (void**)queryDensity.contextPtrs(), density.localVector(), &norm); CPPUNIT_ASSERT(!err);
        queryDensity.closeDB(_auxDB);

        const PylithReal tolerance = 1.0e-6;
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, norm, tolerance);
    } // Test getting density field

    { // Test getting bulkModulus field.
        pylith::topology::Field bulkModulus(*_mesh);
        _mymaterial->getAuxField(&bulkModulus, "bulk_modulus");

        //bulkModulus.view("BULK MODULUS"); // DEBUGGING

        // Check result
        CPPUNIT_ASSERT_EQUAL(std::string("bulk_modulus"), std::string(bulkModulus.label()));
        CPPUNIT_ASSERT_EQUAL(_mydata->dimension, bulkModulus.spaceDim());

        pylith::topology::FieldQuery queryBulkModulus(bulkModulus);
        queryBulkModulus.queryFn("bulk_modulus", pylith::materials::Query::dbQueryBulkModulus);
        queryBulkModulus.openDB(_auxDB, lengthScale);

        PylithReal norm = 0.0;
        const PylithReal t = _mydata->t;
        const PetscDM dm = bulkModulus.dmMesh(); CPPUNIT_ASSERT(dm);
        PetscErrorCode err = DMComputeL2Diff(dm, t, queryBulkModulus.functions(), (void**)queryBulkModulus.contextPtrs(), bulkModulus.localVector(), &norm); CPPUNIT_ASSERT(!err);
        queryBulkModulus.closeDB(_auxDB);

        const PylithReal tolerance = 1.0e-6;
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, norm, tolerance);
    } // Test getting bulkModulus field

    { // Test getting maxwellTime field.
        pylith::topology::Field maxwellTime(*_mesh);
        _mymaterial->getAuxField(&maxwellTime, "maxwell_time");

        //maxwellTime.view("MAXWELL TIME"); // DEBUGGING

        // Check result
        CPPUNIT_ASSERT_EQUAL(std::string("maxwell_time"), std::string(maxwellTime.label()));
        CPPUNIT_ASSERT_EQUAL(_mydata->dimension, maxwellTime.spaceDim());

        pylith::topology::FieldQuery queryMaxwellTime(maxwellTime);
        queryMaxwellTime.queryFn("maxwell_time", pylith::materials::Query::dbQueryMaxwellTime);
        queryMaxwellTime.openDB(_auxDB, lengthScale);

        PylithReal norm = 0.0;
        const PylithReal t = _mydata->t;
        const PetscDM dm = maxwellTime.dmMesh(); CPPUNIT_ASSERT(dm);
        PetscErrorCode err = DMComputeL2Diff(dm, t, queryMaxwellTime.functions(), (void**)queryMaxwellTime.contextPtrs(), maxwellTime.localVector(), &norm); CPPUNIT_ASSERT(!err);
        queryMaxwellTime.closeDB(_auxDB);

        const PylithReal tolerance = 1.0e-6;
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, norm, tolerance);
    } // Test getting maxwellTime field

    std::cout << " :TODO: Add test for getting reference stress/strain. " << std::endl;

    PYLITH_METHOD_END;
} // testGetAuxField


// ----------------------------------------------------------------------
// Get material.
pylith::materials::MaterialNew*
pylith::materials::TestIsotropicLinearMaxwellPlaneStrain::_material(void)
{ // _material
    return _mymaterial;
} // _material


// ----------------------------------------------------------------------
// Get test data.
pylith::materials::TestMaterialNew_Data*
pylith::materials::TestIsotropicLinearMaxwellPlaneStrain::_data(void)
{ // _data
    return _mydata;
} // _data


// ----------------------------------------------------------------------
// Setup and populate solution fields.
void
pylith::materials::TestIsotropicLinearMaxwellPlaneStrain::_setupSolutionFields(void)
{ // _setupSolutionFields
    PYLITH_METHOD_BEGIN;

    CPPUNIT_ASSERT(_solutionFields);

    CPPUNIT_ASSERT( (!_mydata->isExplicit && 1 == _mydata->numSolnFields) ||
                    (_mydata->isExplicit && 2 == _mydata->numSolnFields) );
    CPPUNIT_ASSERT(_mydata->solnDiscretizations);
    CPPUNIT_ASSERT(_mydata->normalizer);

    spatialdata::spatialdb::SimpleGridDB solutionDB;
    solutionDB.filename(_mydata->solnDBFilename);
    solutionDB.label("solution database");
    solutionDB.queryType(spatialdata::spatialdb::SimpleGridDB::LINEAR);

    { // Solution @ t1
        pylith::topology::Field& solution = _solutionFields->get("solution");
        pylith::problems::SolutionFactory factory(solution, *_mydata->normalizer);
        factory.displacement(_mydata->solnDiscretizations[0]);
        if (_mydata->isExplicit) {
            factory.velocity(_mydata->solnDiscretizations[1]);
        } // if
        solution.subfieldsSetup();
        solution.allocate();
        factory.setValues(&solutionDB);
    } // Solution @ t1

    { // Time derivative of solution @ t1
        pylith::topology::Field& solutionDot = _solutionFields->get("solution_dot");
        pylith::problems::SolutionFactory factory(solutionDot, *_mydata->normalizer);
        factory.displacementDot(_mydata->solnDiscretizations[0]);
        if (_mydata->isExplicit) {
            factory.velocityDot(_mydata->solnDiscretizations[1]);
        } // if
        solutionDot.subfieldsSetup();
        solutionDot.allocate();
        factory.setValues(&solutionDB);
    } // Time derivative of solution @ t1

    spatialdata::spatialdb::SimpleGridDB perturbationDB;
    perturbationDB.filename(_mydata->pertDBFilename);
    perturbationDB.label("perturbation database");
    perturbationDB.queryType(spatialdata::spatialdb::SimpleGridDB::LINEAR);

    { // Perturbation @ t2
        pylith::topology::Field& perturbation = _solutionFields->get("perturbation");
        pylith::problems::SolutionFactory factory(perturbation, *_mydata->normalizer);
        perturbation.cloneSection(_solutionFields->get("solution"));
        perturbation.allocate();
        factory.setValues(&perturbationDB);
    } // Perturbation @ t2

    { // Time derivative of solution @ t2
        pylith::topology::Field& perturbationDot = _solutionFields->get("perturbation_dot");
        pylith::problems::SolutionFactory factory(perturbationDot, *_mydata->normalizer);
        perturbationDot.cloneSection(_solutionFields->get("solution_dot"));
        perturbationDot.allocate();
        factory.setValues(&perturbationDB);
    } // Time derivative of solution @ t2

} // _setupSolutionFields

// ----------------------------------------------------------------------
// Constructor
pylith::materials::TestIsotropicLinearMaxwellPlaneStrain_Data::TestIsotropicLinearMaxwellPlaneStrain_Data(void) :
    useInertia(false),
    useBodyForce(false),
    useGravity(false),
    useReferenceState(false)
{ // constructor
    dimension = 2;
    gravityVector[0] = 0.0; // Use scales in test to provide correct nondimensional value.
    gravityVector[1] = 0.0;
    gravityVector[2] = 0;
} // constructor


// ----------------------------------------------------------------------
// Destructor
pylith::materials::TestIsotropicLinearMaxwellPlaneStrain_Data::~TestIsotropicLinearMaxwellPlaneStrain_Data(void)
{ // destructor
} // destructor


// End of file
