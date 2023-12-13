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

#include "TestAuxiliaryFactoryWellboreSource.hh" // Implementation of class methods

#include "pylith/sources/AuxiliaryFactoryWellboreSource.hh" // Test subject
#include "tests/src/FieldTester.hh" // USES FieldTester

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/MeshOps.hh" // USES MeshOps
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/meshio/MeshIOAscii.hh" // USES MeshIOAscii

#include "pylith/utils/error.hh" // USES PYLITH_METHOD*
#include "pylith/utils/journals.hh" // USES pythia::journal::debug_t

#include "spatialdata/spatialdb/UserFunctionDB.hh" // USES UserFunctionDB
#include "spatialdata/spatialdb/GravityField.hh" // USES GravityField
#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include "catch2/catch_test_macros.hpp"

// ------------------------------------------------------------------------------------------------
// Constructor.
pylith::sources::TestAuxiliaryFactoryWellboreSource::TestAuxiliaryFactoryWellboreSource(TestAuxiliaryFactoryWellboreSource_Data* data) :
    _data(data) {
    PYLITH_METHOD_BEGIN;

    assert(_data->normalizer);
    _data->normalizer->setLengthScale(1.0e+03);
    _data->normalizer->setTimeScale(2.0);
    _data->normalizer->setDensityScale(3.0e+3);
    _data->normalizer->setPressureScale(2.25e+10);

    pylith::topology::Field::SubfieldInfo info;
    pylith::string_vector componentNames;

    // fluid_density
    componentNames.resize(1);
    componentNames[0] = "fluid_density";
    info.description = pylith::topology::Field::Description(
        "fluid_density",
        "fluid_density",
        componentNames,
        componentNames.size(),
        pylith::topology::Field::SCALAR,
        _data->normalizer->getDensityScale(),
        pylith::topology::FieldQuery::validatorPositive
        );
    info.fe = pylith::topology::Field::Discretization(
        1, 2, _data->auxDim, 1, false, pylith::topology::Field::DEFAULT_BASIS, pylith::topology::Field::POLYNOMIAL_SPACE, true
        );
    info.index = 0;
    _data->subfields["fluid_density"] = info;

    // fluid_viscosity
    componentNames.resize(1);
    componentNames[0] = "fluid_viscosity";
    info.description = pylith::topology::Field::Description(
        "fluid_viscosity",
        "fluid_viscosity",
        componentNames,
        componentNames.size(),
        pylith::topology::Field::SCALAR,
        _data->normalizer->getViscosityScale(),
        pylith::topology::FieldQuery::validatorPositive
        );
    info.fe = pylith::topology::Field::Discretization(
        1, 2, _data->auxDim, 1, false, pylith::topology::Field::DEFAULT_BASIS, pylith::topology::Field::POLYNOMIAL_SPACE, true
        );
    info.index = 0;
    _data->subfields["fluid_viscosity"] = info;

    // isotropic_permeability
    componentNames.resize(1);
    componentNames[0] = "isotropic_permeability";
    info.description = pylith::topology::Field::Description(
        "isotropic_permeability",
        "isotropic_permeability",
        componentNames,
        componentNames.size(),
        pylith::topology::Field::SCALAR,
        _data->normalizer->getLengthScale() * normalizer->getLengthScale(),
        pylith::topology::FieldQuery::validatorPositive
        );
    info.fe = pylith::topology::Field::Discretization(
        1, 2, _data->auxDim, 1, false, pylith::topology::Field::DEFAULT_BASIS, pylith::topology::Field::POLYNOMIAL_SPACE, true
        );
    info.index = 0;
    _data->subfields["isotropic_permeability"] = info;

    // wellbore_radius
    componentNames.resize(1);
    componentNames[0] = "wellbore_radius";
    info.description = pylith::topology::Field::Description(
        "wellbore_radius",
        "wellbore_radius",
        componentNames,
        componentNames.size(),
        pylith::topology::Field::SCALAR,
        _data->normalizer->getLengthScale(),
        pylith::topology::FieldQuery::validatorPositive
        );
    info.fe = pylith::topology::Field::Discretization(
        1, 2, _data->auxDim, 1, false, pylith::topology::Field::DEFAULT_BASIS, pylith::topology::Field::POLYNOMIAL_SPACE, true
        );
    info.index = 0;
    _data->subfields["wellbore_radius"] = info;

    // wellbore_length
    componentNames.resize(1);
    componentNames[0] = "wellbore_length";
    info.description = pylith::topology::Field::Description(
        "wellbore_length",
        "wellbore_length",
        componentNames,
        componentNames.size(),
        pylith::topology::Field::SCALAR,
        _data->normalizer->getLengthScale(),
        pylith::topology::FieldQuery::validatorPositive
        );
    info.fe = pylith::topology::Field::Discretization(
        1, 2, _data->auxDim, 1, false, pylith::topology::Field::DEFAULT_BASIS, pylith::topology::Field::POLYNOMIAL_SPACE, true
        );
    info.index = 0;
    _data->subfields["wellbore_length"] = info;

    // wellbore_pressure
    componentNames.resize(1);
    componentNames[0] = "wellbore_pressure";
    info.description = pylith::topology::Field::Description(
        "wellbore_pressure",
        "wellbore_pressure",
        componentNames,
        componentNames.size(),
        pylith::topology::Field::SCALAR,
        _data->normalizer->getPressureScale(),
        pylith::topology::FieldQuery::validatorPositive
        );
    info.fe = pylith::topology::Field::Discretization(
        1, 2, _data->auxDim, 1, false, pylith::topology::Field::DEFAULT_BASIS, pylith::topology::Field::POLYNOMIAL_SPACE, true
        );
    info.index = 0;
    _data->subfields["wellbore_pressure"] = info;

    // wellbore_character
    componentNames.resize(1);
    componentNames[0] = "wellbore_character";
    info.description = pylith::topology::Field::Description(
        "wellbore_character",
        "wellbore_character",
        componentNames,
        componentNames.size(),
        pylith::topology::Field::SCALAR,
        1,
        NULL
        );
    info.fe = pylith::topology::Field::Discretization(
        1, 2, _data->auxDim, 1, false, pylith::topology::Field::DEFAULT_BASIS, pylith::topology::Field::POLYNOMIAL_SPACE, true
        );
    info.index = 0;
    _data->subfields["wellbore_character"] = info;

    // element_dimensions
    componentNames.resize(3);
    componentNames[0] = "element_dimensions_x";
    componentNames[1] = "element_dimensions_y";
    componentNames[2] = "element_dimensions_z";

    info.description = pylith::topology::Field::Description(
        "element_dimensions",
        "element_dimensions",
        componentNames,
        componentNames.size(),
        pylith::topology::Field::TENSOR,
        _data->normalizer->getLengthScale(),
        pylith::topology::FieldQuery::validatorPositive
        );
    info.fe = pylith::topology::Field::Discretization(
        2, 2, _data->auxDim, _data->auxDim, false, pylith::topology::Field::DEFAULT_BASIS, pylith::topology::Field::POLYNOMIAL_SPACE, false
        );
    info.index = 0;
    _data->subfields["element_dimensions"] = info;
    if (2 == _data->auxDim) {
        _data->subfields["element_dimensions"].description.numComponents = 2;
        _data->subfields["element_dimensions"].description.vectorFieldType = pylith::topology::Field::OTHER;
    } // if

    // timeDelay
    componentNames.resize(1);
    componentNames[0] = "time_delay";
    info.description = pylith::topology::Field::Description(
        "time_delay",
        "time_delay",
        componentNames,
        componentNames.size(),
        pylith::topology::Field::SCALAR,
        _data->normalizer->getTimeScale(),
        pylith::topology::FieldQuery::validatorNonnegative
        );
    info.fe = pylith::topology::Field::Discretization(
        1, 2, _data->auxDim, 1, false, pylith::topology::Field::DEFAULT_BASIS, pylith::topology::Field::POLYNOMIAL_SPACE, true
        );
    info.index = 1;
    _data->subfields["time_delay"] = info;

    _initialize();

    PYLITH_METHOD_END;
} // constructor


// ------------------------------------------------------------------------------------------------
// Destructor.
pylith::sources::TestAuxiliaryFactoryWellboreSource::~TestAuxiliaryFactoryWellboreSource(void) {
    PYLITH_METHOD_BEGIN;

    delete _factory;_factory = NULL;
    delete _data;_data = NULL;
    delete _mesh;_mesh = NULL;
    delete _auxiliaryField;_auxiliaryField = NULL;

    PYLITH_METHOD_END;
} // tearDown


// ------------------------------------------------------------------------------------------------
// Test adding moment tensor and time delay subfields.
void
pylith::sources::TestAuxiliaryFactoryWellboreSource::testAdd(void) {
    PYLITH_METHOD_BEGIN;
    assert(_factory);
    assert(_data);

    CHECK(!_auxiliaryField->hasSubfield("fluid_density"));
    CHECK(!_auxiliaryField->hasSubfield("fluid_viscosity"));
    CHECK(!_auxiliaryField->hasSubfield("isotropic_permeability"));
    CHECK(!_auxiliaryField->hasSubfield("wellbore_radius"));
    CHECK(!_auxiliaryField->hasSubfield("wellbore_length"));
    CHECK(!_auxiliaryField->hasSubfield("wellbore_pressure"));
    CHECK(!_auxiliaryField->hasSubfield("wellbore_character"));
    CHECK(!_auxiliaryField->hasSubfield("element_dimensions"));
    CHECK(!_auxiliaryField->hasSubfield("time_delay"));

    _factory->addFluidDensity();
    _factory->addFluidViscosity();
    _factory->addIsotropicPermeability();
    _factory->addWellboreRadius();
    _factory->addWellboreLength();
    _factory->addWellboreCharacter();
    _factory->addElementDimensions();
    _factory->addTimeDelay();

    assert(_data->normalizer);

    pylith::testing::FieldTester::checkSubfieldInfo(*_auxiliaryField, _data->subfields["fluid_density"]);
    pylith::testing::FieldTester::checkSubfieldInfo(*_auxiliaryField, _data->subfields["fluid_viscosity"]);
    pylith::testing::FieldTester::checkSubfieldInfo(*_auxiliaryField, _data->subfields["isotropic_permeability"]);
    pylith::testing::FieldTester::checkSubfieldInfo(*_auxiliaryField, _data->subfields["wellbore_radius"]);
    pylith::testing::FieldTester::checkSubfieldInfo(*_auxiliaryField, _data->subfields["wellbore_length"]);
    pylith::testing::FieldTester::checkSubfieldInfo(*_auxiliaryField, _data->subfields["wellbore_pressure"]);
    pylith::testing::FieldTester::checkSubfieldInfo(*_auxiliaryField, _data->subfields["wellbore_character"]);
    pylith::testing::FieldTester::checkSubfieldInfo(*_auxiliaryField, _data->subfields["element_dimensions"]);
    pylith::testing::FieldTester::checkSubfieldInfo(*_auxiliaryField, _data->subfields["time_delay"]);

    PYLITH_METHOD_END;
} // testAdd


// ------------------------------------------------------------------------------------------------
// Test setValues().
void
pylith::sources::TestAuxiliaryFactoryWellboreSource::testSetValuesFromDB(void) {
    PYLITH_METHOD_BEGIN;

    assert(_factory);

    _factory->addFluidDensity();
    _factory->addFluidViscosity();
    _factory->addIsotropicPermeability();
    _factory->addWellboreRadius();
    _factory->addWellboreLength();
    _factory->addWellborePressure();
    _factory->addWellboreCharacter();
    _factory->addElementDimensions();
    _factory->addTimeDelay();
    _auxiliaryField->subfieldsSetup();
    _auxiliaryField->createDiscretization();
    _auxiliaryField->allocate();

    assert(_data);
    assert(_data->normalizer);
    _factory->setValuesFromDB();
    pylith::testing::FieldTester::checkFieldWithDB(*_auxiliaryField, _data->auxiliaryDB, _data->normalizer->getLengthScale());

    PYLITH_METHOD_END;
} // testSetValues


// ------------------------------------------------------------------------------------------------
// Initialze mesh, coordinate system, auxiliary field, and factory.
void
pylith::sources::TestAuxiliaryFactoryWellboreSource::_initialize(void) {
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
    assert(_data->normalizer);
    pylith::topology::MeshOps::nondimensionalize(_mesh, *_data->normalizer);

    _auxiliaryField = new pylith::topology::Field(*_mesh);assert(_auxiliaryField);
    _auxiliaryField->setLabel("auxiliary");

    _factory = new AuxiliaryFactoryWellboreSource();
    assert(_data->auxiliaryDB);
    _factory->setQueryDB(_data->auxiliaryDB);
    typedef std::map<std::string, pylith::topology::Field::SubfieldInfo>::const_iterator subfield_iter;
    for (subfield_iter iter = _data->subfields.begin(); iter != _data->subfields.end(); ++iter) {
        const char* subfieldName = iter->first.c_str();
        const pylith::topology::Field::Discretization& fe = iter->second.fe;
        _factory->setSubfieldDiscretization(subfieldName, fe.basisOrder, fe.quadOrder, fe.dimension, fe.isFaultOnly, fe.cellBasis,
                                            fe.feSpace, fe.isBasisContinuous);
    } // for
    assert(_data->normalizer);
    _factory->initialize(_auxiliaryField, *_data->normalizer, _data->dimension);

    PYLITH_METHOD_END;
} // _initialize


// ------------------------------------------------------------------------------------------------
pylith::sources::TestAuxiliaryFactoryWellboreSource_Data::TestAuxiliaryFactoryWellboreSource_Data(void) :
    meshFilename(NULL),
    cs(NULL),
    normalizer(new spatialdata::units::Nondimensional),
    auxiliaryDB(new spatialdata::spatialdb::UserFunctionDB) {}


// ------------------------------------------------------------------------------------------------
pylith::sources::TestAuxiliaryFactoryWellboreSource_Data::~TestAuxiliaryFactoryWellboreSource_Data(void) {
    delete cs;cs = NULL;
    delete normalizer;normalizer = NULL;
    delete auxiliaryDB;auxiliaryDB = NULL;
} // testAdd


// ------------------------------------------------------------------------------------------------
// Test setValues().
void
pylith::sources::TestAuxiliaryFactoryWellboreSource::testSetValuesFromDB(void) {
    PYLITH_METHOD_BEGIN;

    assert(_factory);

    _factory->addFluidDensity();
    _factory->addFluidViscosity();
    _factory->addIsotropicPermeability();
    _factory->addWellboreRadius();
    _factory->addWellboreLength();
    _factory->addWellborePressure();
    _factory->addWellboreCharacter();
    _factory->addElementDimensions();
    _factory->addTimeDelay();
    _auxiliaryField->subfieldsSetup();
    _auxiliaryField->createDiscretization();
    _auxiliaryField->allocate();

    assert(_data);
    assert(_data->normalizer);
    _factory->setValuesFromDB();
    pylith::testing::FieldTester::checkFieldWithDB(*_auxiliaryField, _data->auxiliaryDB, _data->normalizer->getLengthScale());

    PYLITH_METHOD_END;
} // testSetValues


// ------------------------------------------------------------------------------------------------
// Initialze mesh, coordinate system, auxiliary field, and factory.
void
pylith::sources::TestAuxiliaryFactoryWellboreSource::_initialize(void) {
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
    assert(_data->normalizer);
    pylith::topology::MeshOps::nondimensionalize(_mesh, *_data->normalizer);

    _auxiliaryField = new pylith::topology::Field(*_mesh);assert(_auxiliaryField);
    _auxiliaryField->setLabel("auxiliary");

    _factory = new AuxiliaryFactoryWellboreSource();
    assert(_data->auxiliaryDB);
    _factory->setQueryDB(_data->auxiliaryDB);
    typedef std::map<std::string, pylith::topology::Field::SubfieldInfo>::const_iterator subfield_iter;
    for (subfield_iter iter = _data->subfields.begin(); iter != _data->subfields.end(); ++iter) {
        const char* subfieldName = iter->first.c_str();
        const pylith::topology::Field::Discretization& fe = iter->second.fe;
        _factory->setSubfieldDiscretization(subfieldName, fe.basisOrder, fe.quadOrder, fe.dimension, fe.isFaultOnly, fe.cellBasis,
                                            fe.feSpace, fe.isBasisContinuous);
    } // for
    assert(_data->normalizer);
    _factory->initialize(_auxiliaryField, *_data->normalizer, _data->dimension);

    PYLITH_METHOD_END;
} // _initialize


// ------------------------------------------------------------------------------------------------
pylith::sources::TestAuxiliaryFactoryWellboreSource_Data::TestAuxiliaryFactoryWellboreSource_Data(void) :
    meshFilename(NULL),
    cs(NULL),
    normalizer(new spatialdata::units::Nondimensional),
    auxiliaryDB(new spatialdata::spatialdb::UserFunctionDB) {}


// ------------------------------------------------------------------------------------------------
pylith::sources::TestAuxiliaryFactoryWellboreSource_Data::~TestAuxiliaryFactoryWellboreSource_Data(void) {
    delete cs;cs = NULL;
    delete normalizer;normalizer = NULL;
    delete auxiliaryDB;auxiliaryDB = NULL;
} // destructor
