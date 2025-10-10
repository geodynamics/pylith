// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================

#include <portinfo>

#include "TestAuxiliaryFactoryElasticity.hh" // Implementation of class methods

#include "pylith/materials/AuxiliaryFactoryElasticity.hh" // Test subject
#include "tests/src/FieldTester.hh" // USES FieldTester

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/MeshOps.hh" // USES MeshOps
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/meshio/MeshIOAscii.hh" // USES MeshIOAscii
#include "pylith/utils/error.hh" // USES PYLITH_METHOD_*

#include "spatialdata/spatialdb/UserFunctionDB.hh" // USES UserFunctionDB
#include "spatialdata/spatialdb/GravityField.hh" // USES GravityField
#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys
#include "pylith/scales/Scales.hh" // USES Scales
#include "pylith/scales/ElasticityScales.hh" // USES ElasticityScales

#include "catch2/catch_test_macros.hpp"

// ------------------------------------------------------------------------------------------------
// Constructor.
pylith::materials::TestAuxiliaryFactoryElasticity::TestAuxiliaryFactoryElasticity(TestAuxiliaryFactoryElasticity_Data* data) :
    _data(data) {
    PYLITH_METHOD_BEGIN;

    pylith::topology::Field::SubfieldInfo info;
    pylith::string_vector componentNames;

    // density
    componentNames.resize(1);
    componentNames[0] = "density";
    info.description = pylith::topology::Field::Description(
        "density",
        "density",
        componentNames,
        componentNames.size(),
        pylith::topology::Field::SCALAR,
        pylith::scales::ElasticityScales::getDensityScale(*_data->scales),
        0.0,
        pylith::topology::FieldQuery::validatorPositive
        );
    info.fe = pylith::topology::Field::Discretization(
        1, 2, _data->auxDim, 1, false, pylith::topology::Field::DEFAULT_BASIS, pylith::topology::Field::POLYNOMIAL_SPACE, true
        );
    info.index = 0;
    _data->subfields["density"] = info;

    // body_force
    componentNames.resize(3);
    componentNames[0] = "body_force_x";
    componentNames[1] = "body_force_y";
    componentNames[2] = "body_force_z";
    info.description = pylith::topology::Field::Description(
        "body_force",
        "body_force",
        componentNames,
        componentNames.size(),
        pylith::topology::Field::VECTOR,
        pylith::scales::ElasticityScales::getBodyForceScale(*_data->scales)
        );
    info.fe = pylith::topology::Field::Discretization(
        2, 2, _data->auxDim, _data->auxDim, false, pylith::topology::Field::DEFAULT_BASIS, pylith::topology::Field::POLYNOMIAL_SPACE, false
        );
    info.index = 1;
    _data->subfields["body_force"] = info;
    _data->subfields["body_force"].description.numComponents = _data->auxDim;

    // gravity_field
    componentNames.resize(3);
    componentNames[0] = "gravitational_acceleration_x";
    componentNames[1] = "gravitational_acceleration_y";
    componentNames[2] = "gravitational_acceleration_z";
    info.description = pylith::topology::Field::Description(
        "gravitational_acceleration",
        "gravitational_acceleration",
        componentNames,
        componentNames.size(),
        pylith::topology::Field::VECTOR,
        pylith::scales::ElasticityScales::getAccelerationScale(*_data->scales)
        );
    info.fe = pylith::topology::Field::Discretization(
        2, 2, _data->auxDim, _data->auxDim, false, pylith::topology::Field::DEFAULT_BASIS, pylith::topology::Field::POLYNOMIAL_SPACE, true
        );
    info.index = 2;
    _data->subfields["gravitational_acceleration"] = info;
    _data->subfields["gravitational_acceleration"].description.numComponents = _data->auxDim;

    _initialize();

    PYLITH_METHOD_END;
} // constructor


// ------------------------------------------------------------------------------------------------
// Destructor.
pylith::materials::TestAuxiliaryFactoryElasticity::~TestAuxiliaryFactoryElasticity(void) {
    PYLITH_METHOD_BEGIN;

    delete _factory;_factory = NULL;
    delete _data;_data = NULL;
    delete _mesh;_mesh = NULL;
    delete _auxiliaryField;_auxiliaryField = NULL;

    PYLITH_METHOD_END;
} // tearDown


// ------------------------------------------------------------------------------------------------
// Test adding density, body force, and gravity subfields.
void
pylith::materials::TestAuxiliaryFactoryElasticity::testAdd(void) {
    PYLITH_METHOD_BEGIN;

    assert(_factory);
    assert(_data);

    CHECK(!_auxiliaryField->hasSubfield("density"));
    CHECK(!_auxiliaryField->hasSubfield("body_force"));
    CHECK(!_auxiliaryField->hasSubfield("gravitational_acceleration"));

    _factory->addDensity();
    _factory->addBodyForce();
    assert(_data->gravityField);
    _factory->addGravityField(_data->gravityField);

    assert(_data->scales);

    pylith::testing::FieldTester::checkSubfieldInfo(*_auxiliaryField, _data->subfields["density"]);
    pylith::testing::FieldTester::checkSubfieldInfo(*_auxiliaryField, _data->subfields["body_force"]);
    pylith::testing::FieldTester::checkSubfieldInfo(*_auxiliaryField, _data->subfields["gravitational_acceleration"]);

    PYLITH_METHOD_END;
} // testAdd


// ------------------------------------------------------------------------------------------------
// Test setValues().
void
pylith::materials::TestAuxiliaryFactoryElasticity::testSetValuesFromDB(void) {
    PYLITH_METHOD_BEGIN;

    assert(_factory);

    _factory->addDensity();
    _factory->addBodyForce();
    assert(_data->gravityField);
    _factory->addGravityField(_data->gravityField);
    _auxiliaryField->subfieldsSetup();
    _auxiliaryField->createDiscretization();
    _auxiliaryField->allocate();

    assert(_data);
    assert(_data->scales);
    _factory->setValuesFromDB();
    pylith::testing::FieldTester::checkFieldWithDB(*_auxiliaryField, _data->auxiliaryDB, _data->scales->getLengthScale());

    PYLITH_METHOD_END;
} // testSetValues


// ------------------------------------------------------------------------------------------------
// Initialze mesh, coordinate system, auxiliary field, and factory.
void
pylith::materials::TestAuxiliaryFactoryElasticity::_initialize(void) {
    PYLITH_METHOD_BEGIN;
    assert(_data);

    pylith::meshio::MeshIOAscii iohandler;
    assert(_data->meshFilename);
    iohandler.setFilename(_data->meshFilename);
    _mesh = new pylith::topology::Mesh();assert(_mesh);
    iohandler.read(_mesh);

    assert(pylith::topology::MeshOps::getNumCells(*_mesh) > 0);
    assert(pylith::topology::MeshOps::getNumVertices(*_mesh) > 0);

    // Setup coordinates.
    _mesh->setCoordSys(_data->cs);
    assert(_data->scales);
    pylith::topology::MeshOps::nondimensionalize(_mesh, *_data->scales);

    _auxiliaryField = new pylith::topology::Field(*_mesh);assert(_auxiliaryField);
    _auxiliaryField->setLabel("auxiliary");

    _factory = new AuxiliaryFactoryElasticity();
    assert(_data->auxiliaryDB);
    _factory->setQueryDB(_data->auxiliaryDB);
    typedef std::map<std::string, pylith::topology::Field::SubfieldInfo>::const_iterator subfield_iter;
    for (subfield_iter iter = _data->subfields.begin(); iter != _data->subfields.end(); ++iter) {
        const char* subfieldName = iter->first.c_str();
        const pylith::topology::Field::Discretization& fe = iter->second.fe;
        _factory->setSubfieldDiscretization(subfieldName, fe.basisOrder, fe.quadOrder, fe.dimension, fe.isFaultOnly, fe.cellBasis,
                                            fe.feSpace, fe.isBasisContinuous);
    } // for
    assert(_data->scales);
    _factory->initialize(_auxiliaryField, *_data->scales, _data->dimension);

    PYLITH_METHOD_END;
} // _initialize


// ------------------------------------------------------------------------------------------------
pylith::materials::TestAuxiliaryFactoryElasticity_Data::TestAuxiliaryFactoryElasticity_Data(void) :
    meshFilename(NULL),
    cs(NULL),
    scales(new pylith::scales::Scales),
    auxiliaryDB(new spatialdata::spatialdb::UserFunctionDB),
    gravityField(new spatialdata::spatialdb::GravityField) {}


// ------------------------------------------------------------------------------------------------
pylith::materials::TestAuxiliaryFactoryElasticity_Data::~TestAuxiliaryFactoryElasticity_Data(void) {
    delete cs;cs = NULL;
    delete scales;scales = NULL;
    delete auxiliaryDB;auxiliaryDB = NULL;
    delete gravityField;gravityField = NULL;
} // destructor


// End of file
