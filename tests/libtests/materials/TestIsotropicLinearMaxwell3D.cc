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

#include "TestIsotropicLinearMaxwell3D.hh" // Implementation of class methods

#include "pylith/materials/IsotropicLinearMaxwell3D.hh" // USES IsotropicLinearMaxwell3D
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
pylith::materials::TestIsotropicLinearMaxwell3D::setUp(void) {
    TestMaterial::setUp();
    _mymaterial = new IsotropicLinearMaxwell3D();CPPUNIT_ASSERT(_mymaterial);
    _mydata = NULL;

    GenericComponent::setName("TestIsotropicLinearMaxwell3D");

    _mymaterial->PyreComponent::identifier("TestIsotropicLinearMaxwell3D");
    const char* journal = _mymaterial->PyreComponent::getName();
} // setUp


// ----------------------------------------------------------------------
// Deallocate testing data.
void
pylith::materials::TestIsotropicLinearMaxwell3D::tearDown(void) {
    const char* journal = _mymaterial->PyreComponent::getName();
    pythia::journal::debug_t debug(journal);
    debug.deactivate(); // DEBUGGING

    TestMaterial::tearDown();

    delete _mymaterial;_mymaterial = NULL;
    delete _mydata;_mydata = NULL;
} // tearDown


// ----------------------------------------------------------------------
// Test useInertia(), useBodyForce(), useReferenceState().
void
pylith::materials::TestIsotropicLinearMaxwell3D::testAccessors(void) {
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
pylith::materials::TestIsotropicLinearMaxwell3D::test_auxiliaryFieldSetup(void) {
    PYLITH_METHOD_BEGIN;

    CPPUNIT_ASSERT(_mymaterial);
    CPPUNIT_ASSERT(_mesh);
    CPPUNIT_ASSERT(_mydata);
    CPPUNIT_ASSERT(_mydata->normalizer);
    const PylithReal densityScale = _mydata->normalizer->getDensityScale();
    const PylithReal lengthScale = _mydata->normalizer->getLengthScale();
    const PylithReal timeScale = _mydata->normalizer->getTimeScale();
    const PylithReal pressureScale = _mydata->normalizer->getPressureScale();
    const PylithReal forceScale = pressureScale / lengthScale;
    const PylithReal accelerationScale = lengthScale/(timeScale * timeScale);

    delete _mymaterial->_auxiliaryField;_mymaterial->_auxiliaryField = new topology::Field(*_mesh);CPPUNIT_ASSERT(_mymaterial->_auxiliaryField);
    _mymaterial->_auxiliaryFieldSetup();

    // Check discretizations
    { // density
        const char* label = "density";
        const pylith::topology::Field::SubfieldInfo& info = _mymaterial->_auxiliaryField->subfieldInfo(label);
        CPPUNIT_ASSERT_EQUAL(size_t(1), info.description.numComponents);
        CPPUNIT_ASSERT_EQUAL(std::string(label), info.description.label);
        CPPUNIT_ASSERT_EQUAL(pylith::topology::Field::SCALAR, info.description.vectorFieldType);
        CPPUNIT_ASSERT_EQUAL(densityScale, info.description.scale);
        CPPUNIT_ASSERT_EQUAL(1, info.fe.basisOrder);
        CPPUNIT_ASSERT_EQUAL(1, info.fe.quadOrder);
        CPPUNIT_ASSERT_EQUAL(true, info.fe.isBasisContinuous);
        CPPUNIT_ASSERT_EQUAL(pylith::topology::Field::POLYNOMIAL_SPACE, info.fe.feSpace);
    } // density

    { // shear modulus
        const char* label = "shear_modulus";
        const pylith::topology::Field::SubfieldInfo& info = _mymaterial->_auxiliaryField->subfieldInfo(label);
        CPPUNIT_ASSERT_EQUAL(size_t(1), info.description.numComponents);
        CPPUNIT_ASSERT_EQUAL(std::string(label), info.description.label);
        CPPUNIT_ASSERT_EQUAL(pylith::topology::Field::SCALAR, info.description.vectorFieldType);
        CPPUNIT_ASSERT_EQUAL(pressureScale, info.description.scale);
        CPPUNIT_ASSERT_EQUAL(1, info.fe.basisOrder);
        CPPUNIT_ASSERT_EQUAL(1, info.fe.quadOrder);
        CPPUNIT_ASSERT_EQUAL(true, info.fe.isBasisContinuous);
        CPPUNIT_ASSERT_EQUAL(pylith::topology::Field::POLYNOMIAL_SPACE, info.fe.feSpace);
    } // shear modulus

    { // bulk modulus
        const char* label = "bulk_modulus";
        const pylith::topology::Field::SubfieldInfo& info = _mymaterial->_auxiliaryField->subfieldInfo(label);
        CPPUNIT_ASSERT_EQUAL(size_t(1), info.description.numComponents);
        CPPUNIT_ASSERT_EQUAL(std::string(label), info.description.label);
        CPPUNIT_ASSERT_EQUAL(pylith::topology::Field::SCALAR, info.description.vectorFieldType);
        CPPUNIT_ASSERT_EQUAL(pressureScale, info.description.scale);
        CPPUNIT_ASSERT_EQUAL(1, info.fe.basisOrder);
        CPPUNIT_ASSERT_EQUAL(1, info.fe.quadOrder);
        CPPUNIT_ASSERT_EQUAL(true, info.fe.isBasisContinuous);
        CPPUNIT_ASSERT_EQUAL(pylith::topology::Field::POLYNOMIAL_SPACE, info.fe.feSpace);
    } // bulk modulus

    { // Maxwell time
        const char* label = "maxwell_time";
        const pylith::topology::Field::SubfieldInfo& info = _mymaterial->_auxiliaryField->subfieldInfo(label);
        CPPUNIT_ASSERT_EQUAL(size_t(1), info.description.numComponents);
        CPPUNIT_ASSERT_EQUAL(std::string(label), info.description.label);
        CPPUNIT_ASSERT_EQUAL(pylith::topology::Field::SCALAR, info.description.vectorFieldType);
        CPPUNIT_ASSERT_EQUAL(timeScale, info.description.scale);
        CPPUNIT_ASSERT_EQUAL(1, info.fe.basisOrder);
        CPPUNIT_ASSERT_EQUAL(1, info.fe.quadOrder);
        CPPUNIT_ASSERT_EQUAL(true, info.fe.isBasisContinuous);
        CPPUNIT_ASSERT_EQUAL(pylith::topology::Field::POLYNOMIAL_SPACE, info.fe.feSpace);
    } // Maxwell time

    { // Viscous strain
        const char* label = "viscous_strain";
        const pylith::topology::Field::SubfieldInfo& info = _mymaterial->_auxiliaryField->subfieldInfo(label);
        CPPUNIT_ASSERT_EQUAL(size_t(6), info.description.numComponents);
        CPPUNIT_ASSERT_EQUAL(std::string(label), info.description.label);
        CPPUNIT_ASSERT_EQUAL(pylith::topology::Field::OTHER, info.description.vectorFieldType);
        CPPUNIT_ASSERT_EQUAL(1.0, info.description.scale);
        CPPUNIT_ASSERT_EQUAL(1, info.fe.basisOrder);
        CPPUNIT_ASSERT_EQUAL(1, info.fe.quadOrder);
        CPPUNIT_ASSERT_EQUAL(true, info.fe.isBasisContinuous);
        CPPUNIT_ASSERT_EQUAL(pylith::topology::Field::POLYNOMIAL_SPACE, info.fe.feSpace);
    } // Viscous strain

    { // Total strain
        const char* label = "total_strain";
        const pylith::topology::Field::SubfieldInfo& info = _mymaterial->_auxiliaryField->subfieldInfo(label);
        CPPUNIT_ASSERT_EQUAL(size_t(6), info.description.numComponents);
        CPPUNIT_ASSERT_EQUAL(std::string(label), info.description.label);
        CPPUNIT_ASSERT_EQUAL(pylith::topology::Field::OTHER, info.description.vectorFieldType);
        CPPUNIT_ASSERT_EQUAL(1.0, info.description.scale);
        CPPUNIT_ASSERT_EQUAL(1, info.fe.basisOrder);
        CPPUNIT_ASSERT_EQUAL(1, info.fe.quadOrder);
        CPPUNIT_ASSERT_EQUAL(true, info.fe.isBasisContinuous);
        CPPUNIT_ASSERT_EQUAL(pylith::topology::Field::POLYNOMIAL_SPACE, info.fe.feSpace);
    } // Total strain

    if (_mymaterial->_gravityField) { // gravity field
        const char* label = "gravity_field";
        const pylith::topology::Field::SubfieldInfo& info = _mymaterial->_auxiliaryField->subfieldInfo(label);
        CPPUNIT_ASSERT_EQUAL(size_t(3), info.description.numComponents);
        CPPUNIT_ASSERT_EQUAL(std::string(label), info.description.label);
        CPPUNIT_ASSERT_EQUAL(pylith::topology::Field::VECTOR, info.description.vectorFieldType);
        CPPUNIT_ASSERT_EQUAL(accelerationScale, info.description.scale);
        CPPUNIT_ASSERT_EQUAL(1, info.fe.basisOrder);
        CPPUNIT_ASSERT_EQUAL(1, info.fe.quadOrder);
        CPPUNIT_ASSERT_EQUAL(true, info.fe.isBasisContinuous);
        CPPUNIT_ASSERT_EQUAL(pylith::topology::Field::POLYNOMIAL_SPACE, info.fe.feSpace);
    } // gravity field

    if (_mymaterial->_useBodyForce) { // body force
        const char* label = "body_force";
        const pylith::topology::Field::SubfieldInfo& info = _mymaterial->_auxiliaryField->subfieldInfo(label);
        CPPUNIT_ASSERT_EQUAL(size_t(3), info.description.numComponents);
        CPPUNIT_ASSERT_EQUAL(std::string(label), info.description.label);
        CPPUNIT_ASSERT_EQUAL(pylith::topology::Field::VECTOR, info.description.vectorFieldType);
        CPPUNIT_ASSERT_EQUAL(forceScale, info.description.scale);
        CPPUNIT_ASSERT_EQUAL(1, info.fe.basisOrder);
        CPPUNIT_ASSERT_EQUAL(1, info.fe.quadOrder);
        CPPUNIT_ASSERT_EQUAL(true, info.fe.isBasisContinuous);
        CPPUNIT_ASSERT_EQUAL(pylith::topology::Field::POLYNOMIAL_SPACE, info.fe.feSpace);
    } // body force

    if (_mymaterial->_useReferenceState) { // reference stress
        const char* label = "reference_stress";
        const pylith::topology::Field::SubfieldInfo& info = _mymaterial->_auxiliaryField->subfieldInfo(label);
        CPPUNIT_ASSERT_EQUAL(size_t(6), info.description.numComponents);
        CPPUNIT_ASSERT_EQUAL(std::string(label), info.description.label);
        CPPUNIT_ASSERT_EQUAL(pylith::topology::Field::OTHER, info.description.vectorFieldType);
        CPPUNIT_ASSERT_EQUAL(pressureScale, info.description.scale);
        CPPUNIT_ASSERT_EQUAL(1, info.fe.basisOrder);
        CPPUNIT_ASSERT_EQUAL(1, info.fe.quadOrder);
        CPPUNIT_ASSERT_EQUAL(true, info.fe.isBasisContinuous);
        CPPUNIT_ASSERT_EQUAL(pylith::topology::Field::POLYNOMIAL_SPACE, info.fe.feSpace);
    } // reference stress

    if (_mymaterial->_useReferenceState) { // reference strain
        const char* label = "reference_strain";
        const pylith::topology::Field::SubfieldInfo& info = _mymaterial->_auxiliaryField->subfieldInfo(label);
        CPPUNIT_ASSERT_EQUAL(size_t(6), info.description.numComponents);
        CPPUNIT_ASSERT_EQUAL(std::string(label), info.description.label);
        CPPUNIT_ASSERT_EQUAL(pylith::topology::Field::OTHER, info.description.vectorFieldType);
        CPPUNIT_ASSERT_EQUAL(1.0, info.description.scale);
        CPPUNIT_ASSERT_EQUAL(1, info.fe.basisOrder);
        CPPUNIT_ASSERT_EQUAL(1, info.fe.quadOrder);
        CPPUNIT_ASSERT_EQUAL(true, info.fe.isBasisContinuous);
        CPPUNIT_ASSERT_EQUAL(pylith::topology::Field::POLYNOMIAL_SPACE, info.fe.feSpace);
    } // reference strain

    PYLITH_METHOD_END;
} // test_auxiliaryFieldSetup


// ----------------------------------------------------------------------
// Test getAuxField().
void
pylith::materials::TestIsotropicLinearMaxwell3D::testGetAuxField(void) {
    PYLITH_METHOD_BEGIN;

    _initializeFull();

    CPPUNIT_ASSERT(_mymaterial);
    CPPUNIT_ASSERT(_mesh);
    CPPUNIT_ASSERT(_mydata->normalizer);
    const PylithReal lengthScale = _mydata->normalizer->getLengthScale();

    const pylith::topology::Field* auxField = _mymaterial->auxField();CPPUNIT_ASSERT(auxField);
    { // Test getting density field.
        pylith::topology::Field density(*_mesh);
        density.copySubfield(*auxField, "density");

        // density.view("DENSITY"); // DEBUGGING

        // Check result
        CPPUNIT_ASSERT_EQUAL(std::string("density"), std::string(density.getLabel()));
        CPPUNIT_ASSERT_EQUAL(_mydata->dimension, density.getSpaceDim());

        pylith::topology::FieldQuery queryDensity(density);
        queryDensity.initializeWithDefaultQueryFns();
        queryDensity.openDB(_mydata->auxDB, lengthScale);

        PylithReal norm = 0.0;
        const PylithReal t = _mydata->t;
        const PetscDM dm = density.dmMesh();CPPUNIT_ASSERT(dm);
        PetscErrorCode err = DMPlexComputeL2DiffLocal(dm, t, queryDensity.functions(), (void**)queryDensity.contextPtrs(), density.localVector(), &norm);CPPUNIT_ASSERT(!err);
        queryDensity.closeDB(_mydata->auxDB);

        const PylithReal tolerance = 1.0e-6;
        CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("Test extracting density subfield from auxiliary field failed.", 0.0, norm, tolerance);
    } // Test getting density field

    { // Test getting bulk_modulus field.
        pylith::topology::Field bulkModulus(*_mesh);
        bulkModulus.copySubfield(*auxField, "bulk_modulus");

        // bulkModulus.view("BULK MODULUS"); // DEBUGGING

        // Check result
        CPPUNIT_ASSERT_EQUAL(std::string("bulk_modulus"), std::string(bulkModulus.getLabel()));
        CPPUNIT_ASSERT_EQUAL(_mydata->dimension, bulkModulus.getSpaceDim());

        pylith::topology::FieldQuery queryBulkModulus(bulkModulus);
        queryBulkModulus.initializeWithDefaultQueryFns();
        queryBulkModulus.openDB(_mydata->auxDB, lengthScale);

        PylithReal norm = 0.0;
        const PylithReal t = _mydata->t;
        const PetscDM dm = bulkModulus.dmMesh();CPPUNIT_ASSERT(dm);
        PetscErrorCode err = DMPlexComputeL2DiffLocal(dm, t, queryBulkModulus.functions(), (void**)queryBulkModulus.contextPtrs(), bulkModulus.localVector(), &norm);CPPUNIT_ASSERT(!err);
        queryBulkModulus.closeDB(_mydata->auxDB);

        const PylithReal tolerance = 1.0e-6;
        CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("Test extracting bulk modulus subfield from auxiliary field failed.", 0.0, norm, tolerance);
    } // Test getting bulk_modulus field

    { // Test getting maxwell_time field.
        pylith::topology::Field maxwellTime(*_mesh);
        maxwellTime.copySubfield(*auxField, "maxwell_time");

        // maxwellTime.view("MAXWELL TIME"); // DEBUGGING

        // Check result
        CPPUNIT_ASSERT_EQUAL(std::string("maxwell_time"), std::string(maxwellTime.getLabel()));
        CPPUNIT_ASSERT_EQUAL(_mydata->dimension, maxwellTime.getSpaceDim());

        pylith::topology::FieldQuery queryMaxwellTime(maxwellTime);
        queryMaxwellTime.initializeWithDefaultQueryFns();
        queryMaxwellTime.openDB(_mydata->auxDB, lengthScale);

        PylithReal norm = 0.0;
        const PylithReal t = _mydata->t;
        const PetscDM dm = maxwellTime.dmMesh();CPPUNIT_ASSERT(dm);
        PetscErrorCode err = DMPlexComputeL2DiffLocal(dm, t, queryMaxwellTime.functions(), (void**)queryMaxwellTime.contextPtrs(), maxwellTime.localVector(), &norm);CPPUNIT_ASSERT(!err);
        queryMaxwellTime.closeDB(_mydata->auxDB);

        const PylithReal tolerance = 1.0e-6;
        CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("Test extracting Maxwell time subfield from auxiliary field failed.", 0.0, norm, tolerance);
    } // Test getting maxwell_time field

    { // Test getting viscous_strain field.
        pylith::topology::Field viscousStrain(*_mesh);
        viscousStrain.copySubfield(*auxField, "viscous_strain");

        // viscousStrain.view("VISCOUS STRAIN"); // DEBUGGING

        // Check result
        CPPUNIT_ASSERT_EQUAL(std::string("viscous_strain"), std::string(viscousStrain.getLabel()));
        CPPUNIT_ASSERT_EQUAL(_mydata->dimension, viscousStrain.getSpaceDim());

        pylith::topology::FieldQuery queryViscousStrain(viscousStrain);
        queryViscousStrain.initializeWithDefaultQueryFns();
        queryViscousStrain.openDB(_mydata->auxDB, lengthScale);

        PylithReal norm = 0.0;
        const PylithReal t = _mydata->t;
        const PetscDM dm = viscousStrain.dmMesh();CPPUNIT_ASSERT(dm);
        PetscErrorCode err = DMPlexComputeL2DiffLocal(dm, t, queryViscousStrain.functions(), (void**)queryViscousStrain.contextPtrs(), viscousStrain.localVector(), &norm);CPPUNIT_ASSERT(!err);
        queryViscousStrain.closeDB(_mydata->auxDB);

        const PylithReal tolerance = 1.0e-6;
        CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("Test extracting viscous strain subfield from auxiliary field failed.", 0.0, norm, tolerance);
    } // Test getting viscous_strain field

    { // Test getting total_strain field.
        pylith::topology::Field totalStrain(*_mesh);
        totalStrain.copySubfield(*auxField, "total_strain");

        // totalStrain.view("TOTAL STRAIN"); // DEBUGGING

        // Check result
        CPPUNIT_ASSERT_EQUAL(std::string("total_strain"), std::string(totalStrain.getLabel()));
        CPPUNIT_ASSERT_EQUAL(_mydata->dimension, totalStrain.getSpaceDim());

        pylith::topology::FieldQuery queryTotalStrain(totalStrain);
        queryTotalStrain.initializeWithDefaultQueryFns();
        queryTotalStrain.openDB(_mydata->auxDB, lengthScale);

        PylithReal norm = 0.0;
        const PylithReal t = _mydata->t;
        const PetscDM dm = totalStrain.dmMesh();CPPUNIT_ASSERT(dm);
        PetscErrorCode err = DMPlexComputeL2DiffLocal(dm, t, queryTotalStrain.functions(), (void**)queryTotalStrain.contextPtrs(), totalStrain.localVector(), &norm);CPPUNIT_ASSERT(!err);
        queryTotalStrain.closeDB(_mydata->auxDB);

        const PylithReal tolerance = 1.0e-6;
        CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("Test extracting total strain subfield from auxiliary field failed.", 0.0, norm, tolerance);
    } // Test getting total_strain field

    if (_mymaterial->_useReferenceState) { // Test getting reference_stress field.
        pylith::topology::Field referenceStress(*_mesh);
        referenceStress.copySubfield(*auxField, "reference_stress");

        // referenceStress.view("REFERENCE STRESS"); // DEBUGGING

        // Check result
        CPPUNIT_ASSERT_EQUAL(std::string("reference_stress"), std::string(referenceStress.getLabel()));
        CPPUNIT_ASSERT_EQUAL(_mydata->dimension, referenceStress.getSpaceDim());

        pylith::topology::FieldQuery queryRefStress(referenceStress);
        queryRefStress.initializeWithDefaultQueryFns();
        queryRefStress.openDB(_mydata->auxDB, lengthScale);

        PylithReal norm = 0.0;
        const PylithReal t = _mydata->t;
        const PetscDM dm = referenceStress.dmMesh();CPPUNIT_ASSERT(dm);
        PetscErrorCode err = DMPlexComputeL2DiffLocal(dm, t, queryRefStress.functions(), (void**)queryRefStress.contextPtrs(), referenceStress.localVector(), &norm);CPPUNIT_ASSERT(!err);
        queryRefStress.closeDB(_mydata->auxDB);

        const PylithReal tolerance = 1.0e-6;
        CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("Test extracting reference stress subfield from auxiliary field failed.", 0.0, norm, tolerance);
    } // Test getting reference_stress field

    if (_mymaterial->_useReferenceState) { // Test getting reference_strain field.
        pylith::topology::Field referenceStrain(*_mesh);
        referenceStrain.copySubfield(*auxField, "reference_strain");

        // referenceStrain.view("REFERENCE STRAIN"); // DEBUGGING

        // Check result
        CPPUNIT_ASSERT_EQUAL(std::string("reference_strain"), std::string(referenceStrain.getLabel()));
        CPPUNIT_ASSERT_EQUAL(_mydata->dimension, referenceStrain.getSpaceDim());

        pylith::topology::FieldQuery queryRefStrain(referenceStrain);
        queryRefStrain.initializeWithDefaultQueryFns();
        queryRefStrain.openDB(_mydata->auxDB, lengthScale);

        PylithReal norm = 0.0;
        const PylithReal t = _mydata->t;
        const PetscDM dm = referenceStrain.dmMesh();CPPUNIT_ASSERT(dm);
        PetscErrorCode err = DMPlexComputeL2DiffLocal(dm, t, queryRefStrain.functions(), (void**)queryRefStrain.contextPtrs(), referenceStrain.localVector(), &norm);CPPUNIT_ASSERT(!err);
        queryRefStrain.closeDB(_mydata->auxDB);

        const PylithReal tolerance = 1.0e-6;
        CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("Test extracting reference strain subfield from auxiliary field failed.", 0.0, norm, tolerance);
    } // Test getting reference_strain field

    PYLITH_METHOD_END;
} // testGetAuxField


// ----------------------------------------------------------------------
// Get material.
pylith::materials::Material*
pylith::materials::TestIsotropicLinearMaxwell3D::_material(void) {
    return _mymaterial;
} // _material


// ----------------------------------------------------------------------
// Get test data.
pylith::materials::TestMaterial_Data*
pylith::materials::TestIsotropicLinearMaxwell3D::_data(void) {
    return _mydata;
} // _data


// ----------------------------------------------------------------------
// Setup and populate solution fields.
void
pylith::materials::TestIsotropicLinearMaxwell3D::_setupSolutionFields(void) {
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
        solution.createDiscretization();
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
        solutionDot.createDiscretization();
        solutionDot.allocate();
        factory.setValues(_mydata->solnDB);
    } // Time derivative of solution

    { // Perturbation
        pylith::topology::Field& perturbation = _solutionFields->get("perturbation");
        const pylith::topology::Field& solution = _solutionFields->get("solution");
        perturbation.cloneSection(solution);
        preturbation.createDiscretization();
        perturbation.allocate();
        perturbation.zeroLocal();
        pylith::problems::SolutionFactory factory(perturbation, *_mydata->normalizer);
        factory.setValues(_mydata->perturbDB);
    } // Perturbation

    { // Time derivative perturbation
        pylith::topology::Field& perturbationDot = _solutionFields->get("perturbation_dot");
        const pylith::topology::Field& solutionDot = _solutionFields->get("solution_dot");
        perturbationDot.cloneSection(solutionDot);
        perturbationDot.createDiscretization();
        perturbationDot.allocate();
        perturbationDot.zeroLocal();
        pylith::problems::SolutionFactory factory(perturbationDot, *_mydata->normalizer);
        factory.setValues(_mydata->perturbDB);
    } // Time derivative perturbation

    PYLITH_METHOD_END;
} // _setupSolutionFields


// ----------------------------------------------------------------------
// Constructor
pylith::materials::TestIsotropicLinearMaxwell3D_Data::TestIsotropicLinearMaxwell3D_Data(void) {
    dimension = 3;
    gravityVector[0] = 0.0; // Use scales in test to provide correct nondimensional value.
    gravityVector[1] = 0.0;
    gravityVector[2] = 0;

    cs = new spatialdata::geocoords::CSCart;CPPUNIT_ASSERT(cs);
    cs->setSpaceDim(dimension);

    // Some auxiliary subfields get updated in updateStateVars().
    auxUpdateDB = new spatialdata::spatialdb::UserFunctionDB;CPPUNIT_ASSERT(auxUpdateDB);

    solnDB->setCoordSys(*cs);
    perturbDB->setCoordSys(*cs);
    auxDB->setCoordSys(*cs);
    auxUpdateDB->setCoordSys(*cs);
} // constructor


// ----------------------------------------------------------------------
// Destructor
pylith::materials::TestIsotropicLinearMaxwell3D_Data::~TestIsotropicLinearMaxwell3D_Data(void) {}


// End of file
