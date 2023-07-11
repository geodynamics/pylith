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
// Copyright (c) 2010-2022 University of California, Davis
//
// See LICENSE.md for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "pylith/utils/GenericComponent.hh" // ISA GenericComponent

#include "pylith/feassemble/AuxiliaryFactory.hh" // Test subject

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/MeshOps.hh" // USES MeshOps
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/meshio/MeshIOAscii.hh" // USES MeshIOAscii
#include "pylith/utils/journals.hh"

#include "spatialdata/spatialdb/UserFunctionDB.hh" // USES UserFunctionDB
#include "spatialdata/geocoords/CSCart.hh" // USES CSCart
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include "catch2/catch_test_macros.hpp"
#include "catch2/matchers/catch_matchers_floating_point.hpp"
#include "catch2/matchers/catch_matchers_exception.hpp"

#include <stdexcept>

namespace pylith {
    namespace feassemble {
        class TestAuxiliaryFactory;
    }
}

// ---------------------------------------------------------------------------------------------------------------------
class pylith::feassemble::TestAuxiliaryFactory : public pylith::utils::GenericComponent {
public:

    /// Setup testing data.
    TestAuxiliaryFactory(void);

    /// Tear down testing data.
    ~TestAuxiliaryFactory(void);

    /// Test setQueryDB() and getQueryDB().
    void testQueryDB(void);

    /// Test setSubfieldDiscretization() and getSubfieldDiscretization().
    void testSubfieldDiscretization(void);

    /// Test initialize().
    void testInitialize(void);

    /// Test setValuesFromDB().
    void testSetValuesFromDB(void);

private:

    pylith::feassemble::AuxiliaryFactory* _factory; ///< Test subject.

public:

    static double density(const double x,
                          const double y) {
        return x+y;
    } // density

    static const char* density_units(void) {
        return "kg/m**3";
    } // density_units

    static double velocity_x(const double x,
                             const double y) {
        return 2.0*x - 0.5*y;
    } // velocity_x

    static double velocity_y(const double x,
                             const double y) {
        return 0.4*x + 0.1*y;
    } // velocity_x

    static const char* velocity_units(void) {
        return "m/s";
    } // velocity_units

}; // class TestAuxiliaryFactory

// ---------------------------------------------------------------------------------------------------------------------
// Setup testing data.
pylith::feassemble::TestAuxiliaryFactory::TestAuxiliaryFactory(void) {
    _factory = new AuxiliaryFactory();
} // constructor


// ---------------------------------------------------------------------------------------------------------------------
// Tear down testing data.
pylith::feassemble::TestAuxiliaryFactory::~TestAuxiliaryFactory(void) {
    delete _factory;_factory = NULL;
} // destructor


// ---------------------------------------------------------------------------------------------------------------------
// Test setQueryDB() and getQueryDB().
void
pylith::feassemble::TestAuxiliaryFactory::testQueryDB(void) {
    spatialdata::spatialdb::UserFunctionDB db;
    const std::string& dbLabel = "test database";
    db.setDescription(dbLabel.c_str());

    assert(_factory);
    assert(!_factory->getQueryDB());
    _factory->setQueryDB(&db);

    const spatialdata::spatialdb::SpatialDB* dbTest = _factory->getQueryDB();
    assert(dbTest);
    REQUIRE(dbLabel == std::string(dbTest->getDescription()));
} // testQueryDB


// ---------------------------------------------------------------------------------------------------------------------
// Test setSubfieldDiscretization() and getSubfieldDiscretization().
void
pylith::feassemble::TestAuxiliaryFactory::testSubfieldDiscretization(void) {
    pylith::topology::FieldBase::Discretization feDisp(2, 2, -1, 2, false, pylith::topology::FieldBase::SIMPLEX_BASIS,
                                                       pylith::topology::FieldBase::POLYNOMIAL_SPACE, true);
    pylith::topology::FieldBase::Discretization feVel(3, 2, 1, 2, true, pylith::topology::FieldBase::SIMPLEX_BASIS,
                                                      pylith::topology::FieldBase::POINT_SPACE, false);

    assert(_factory);
    _factory->setSubfieldDiscretization("displacement", feDisp.basisOrder, feDisp.quadOrder, feDisp.dimension,
                                        feDisp.isFaultOnly, feDisp.cellBasis, feDisp.feSpace, feDisp.isBasisContinuous);
    _factory->setSubfieldDiscretization("velocity", feVel.basisOrder, feVel.quadOrder, feVel.dimension,
                                        feVel.isFaultOnly, feVel.cellBasis, feVel.feSpace, feVel.isBasisContinuous);

    { // Check displacement discretization
        const pylith::topology::FieldBase::Discretization& feTest = _factory->getSubfieldDiscretization("displacement");
        CHECK(feDisp.basisOrder == feTest.basisOrder);
        CHECK(feDisp.quadOrder == feTest.quadOrder);
        CHECK(feDisp.dimension == feTest.dimension);
        CHECK(feDisp.isFaultOnly == feTest.isFaultOnly);
        CHECK(feDisp.cellBasis == feTest.cellBasis);
        CHECK(feDisp.feSpace == feTest.feSpace);
        CHECK(feDisp.isBasisContinuous == feTest.isBasisContinuous);
    } // Check displacement discretization

    { // Check velocity discretization
        const pylith::topology::FieldBase::Discretization& feTest = _factory->getSubfieldDiscretization("velocity");
        CHECK(feVel.basisOrder == feTest.basisOrder);
        CHECK(feVel.quadOrder == feTest.quadOrder);
        CHECK(feVel.dimension == feTest.dimension);
        CHECK(feVel.isFaultOnly == feTest.isFaultOnly);
        CHECK(feVel.cellBasis == feTest.cellBasis);
        CHECK(feVel.feSpace == feTest.feSpace);
        CHECK(feVel.isBasisContinuous == feTest.isBasisContinuous);
    } // Check velocity discretization

    { // default for unknown discretization
        const pylith::topology::FieldBase::Discretization& feTest = _factory->getSubfieldDiscretization("xyz");
        CHECK(1 == feTest.basisOrder);
        CHECK(1 == feTest.quadOrder);
        CHECK(-1 == feTest.dimension);
        CHECK(false == feTest.isFaultOnly);
        CHECK(pylith::topology::FieldBase::DEFAULT_BASIS == feTest.cellBasis);
        CHECK(pylith::topology::FieldBase::POLYNOMIAL_SPACE == feTest.feSpace);
        CHECK(true == feTest.isBasisContinuous);
    }
} // testSubfieldDiscretization


// ---------------------------------------------------------------------------------------------------------------------
// Test initialize().
void
pylith::feassemble::TestAuxiliaryFactory::testInitialize(void) {
    pylith::topology::Mesh mesh;
    pylith::topology::Field field(mesh);

    const int spaceDim = 2;
    spatialdata::units::Nondimensional normalizer;
    normalizer.setLengthScale(10.0);

    pylith::topology::Field::Description description;
    description.label = "displacement";
    description.alias = "disp";
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = "displacement_x";
    description.scale = normalizer.getLengthScale();

    assert(_factory);
    _factory->initialize(&field, normalizer, spaceDim, &description);

    CHECK(&field == _factory->_field);
    CHECK(spaceDim == _factory->_spaceDim);
    CHECK(normalizer.getLengthScale() == _factory->_normalizer->getLengthScale());

    const pylith::topology::Field::Description* descriptionTest = _factory->_defaultDescription;
    CHECK(description.label == descriptionTest->label);
    CHECK(description.alias == descriptionTest->alias);
    CHECK(description.vectorFieldType == descriptionTest->vectorFieldType);
    CHECK(description.numComponents == descriptionTest->numComponents);
    CHECK(description.numComponents == descriptionTest->componentNames.size());
    CHECK(description.componentNames[0] == descriptionTest->componentNames[0]);
    CHECK(description.scale == descriptionTest->scale);
    CHECK(description.validator == descriptionTest->validator);

    assert(_factory->_fieldQuery);
} // testInitialize


// ---------------------------------------------------------------------------------------------------------------------
// Test setValuesFromDB().
void
pylith::feassemble::TestAuxiliaryFactory::testSetValuesFromDB(void) {
    const int spaceDim = 2;
    spatialdata::units::Nondimensional normalizer;
    normalizer.setLengthScale(10.0);
    normalizer.setDensityScale(2.0);

    spatialdata::geocoords::CSCart cs;
    cs.setSpaceDim(spaceDim);

    const int numSubfields = 2;
    static const char* subfields[2] = { "density", "velocity" };
    static const pylith::topology::Field::Discretization subfieldDiscretizations[2] = {
        pylith::topology::Field::Discretization(1, 2), // density
        pylith::topology::Field::Discretization(2, 2), // velocity
    };

    pylith::topology::Field::Description descriptionDensity;
    descriptionDensity.label = "density";
    descriptionDensity.alias = "density";
    descriptionDensity.vectorFieldType = pylith::topology::Field::SCALAR;
    descriptionDensity.numComponents = 1;
    descriptionDensity.componentNames.resize(1);
    descriptionDensity.componentNames[0] = "density";
    descriptionDensity.scale = normalizer.getDensityScale();

    pylith::topology::Field::Description descriptionVelocity;
    descriptionVelocity.label = "velocity";
    descriptionVelocity.alias = "velocity";
    descriptionVelocity.vectorFieldType = pylith::topology::Field::VECTOR;
    descriptionVelocity.numComponents = 2;
    descriptionVelocity.componentNames.resize(2);
    descriptionVelocity.componentNames[0] = "velocity_x";
    descriptionVelocity.componentNames[1] = "velocity_y";
    descriptionVelocity.scale = normalizer.getLengthScale() / normalizer.getTimeScale();

    pylith::topology::Field::Description subfieldDescriptions[2];
    subfieldDescriptions[0] = descriptionDensity;
    subfieldDescriptions[1] = descriptionVelocity;

    spatialdata::spatialdb::UserFunctionDB auxiliaryDB;
    auxiliaryDB.addValue("density", TestAuxiliaryFactory::density, TestAuxiliaryFactory::density_units());
    auxiliaryDB.addValue("velocity_x", TestAuxiliaryFactory::velocity_x, TestAuxiliaryFactory::velocity_units());
    auxiliaryDB.addValue("velocity_y", TestAuxiliaryFactory::velocity_y, TestAuxiliaryFactory::velocity_units());
    auxiliaryDB.setCoordSys(cs);

    pylith::topology::Mesh mesh;
    pylith::meshio::MeshIOAscii iohandler;
    iohandler.setFilename("data/tri.mesh");
    iohandler.read(&mesh);
    mesh.setCoordSys(&cs);
    pylith::topology::MeshOps::nondimensionalize(&mesh, normalizer);

    assert(pylith::topology::MeshOps::getNumCells(mesh) > 0);
    assert(pylith::topology::MeshOps::getNumVertices(mesh) > 0);

    assert(_factory);
    _factory->setQueryDB(&auxiliaryDB);

    pylith::topology::Field auxiliaryField(mesh);
    _factory->initialize(&auxiliaryField, normalizer, spaceDim);
    for (int i = 0; i < numSubfields; ++i) {
        auxiliaryField.subfieldAdd(subfieldDescriptions[i], subfieldDiscretizations[i]);
        _factory->setSubfieldQuery(subfields[i]);
    } // for
    auxiliaryField.subfieldsSetup();
    auxiliaryField.createDiscretization();
    auxiliaryField.allocate();

    _factory->setValuesFromDB();

    // Verify auxiliary field
    PylithReal norm = 0.0;
    PylithReal t = 0.0;
    const PetscDM dmField = auxiliaryField.getDM();assert(dmField);
    pylith::topology::FieldQuery query(auxiliaryField);
    query.initializeWithDefaultQueries();
    query.openDB(&auxiliaryDB, normalizer.getLengthScale());
    PetscErrorCode err = DMPlexComputeL2DiffLocal(dmField, t, query._functions, (void**)query._contextPtrs,
                                                  auxiliaryField.getLocalVector(), &norm);assert(!err);
    query.closeDB(&auxiliaryDB);
    const PylithReal tolerance = 1.0e-6;
    CHECK_THAT(norm, Catch::Matchers::WithinAbs(0.0, tolerance));

    AuxiliaryFactory emptyFactory;
    pythia::journal::error_t error("auxiliaryfactory");
    error.deactivate();
    REQUIRE_THROWS_AS(emptyFactory.setValuesFromDB(), std::logic_error);
} // testSetValuesFromDB


// ------------------------------------------------------------------------------------------------
TEST_CASE("TestAuxiliaryFactory::testQueryDB", "[TestAuxiliaryFactory]") {
    pylith::feassemble::TestAuxiliaryFactory().testQueryDB();
}
TEST_CASE("TestAuxiliaryFactory::testSubfieldDiscretization", "[TestAuxiliaryFactory]") {
    pylith::feassemble::TestAuxiliaryFactory().testSubfieldDiscretization();
}
TEST_CASE("TestAuxiliaryFactory::testInitialize", "[TestAuxiliaryFactory]") {
    pylith::feassemble::TestAuxiliaryFactory().testInitialize();
}
TEST_CASE("TestAuxiliaryFactory::testSetValuesFromDB", "[TestAuxiliaryFactory]") {
    pylith::feassemble::TestAuxiliaryFactory().testSetValuesFromDB();
}

// End of file
