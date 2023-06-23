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

#include "TestAuxiliaryFactoryLinearElastic.hh" // Implementation of class methods

#include "pylith/materials/AuxiliaryFactoryElastic.hh" // Test subject
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
pylith::materials::TestAuxiliaryFactoryLinearElastic::TestAuxiliaryFactoryLinearElastic(TestAuxiliaryFactoryLinearElastic_Data* data) :
    _data(data) {
    PYLITH_METHOD_BEGIN;

    assert(_data->normalizer);
    _data->normalizer->setLengthScale(1.0e+03);
    _data->normalizer->setTimeScale(2.0);
    _data->normalizer->setDensityScale(3.0e+3);
    _data->normalizer->setPressureScale(2.25e+10);

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
        _data->normalizer->getDensityScale(),
        pylith::topology::FieldQuery::validatorPositive
        );
    info.fe = pylith::topology::Field::Discretization(
        2, 2, _data->auxDim, 1, false, pylith::topology::Field::DEFAULT_BASIS, pylith::topology::Field::POLYNOMIAL_SPACE, true
        );
    info.index = 0;
    _data->subfields["density"] = info;

    // shear_modulus
    componentNames.resize(1);
    componentNames[0] = "shear_modulus";
    info.description = pylith::topology::Field::Description(
        "shear_modulus",
        "shear_modulus",
        componentNames,
        componentNames.size(),
        pylith::topology::Field::SCALAR,
        _data->normalizer->getPressureScale(),
        pylith::topology::FieldQuery::validatorNonnegative
        );
    info.fe = pylith::topology::Field::Discretization(
        1, 2, _data->auxDim, 1, false, pylith::topology::Field::DEFAULT_BASIS, pylith::topology::Field::POLYNOMIAL_SPACE, true
        );
    info.index = 1;
    _data->subfields["shear_modulus"] = info;

    // bulk_modulus
    componentNames.resize(1);
    componentNames[0] = "bulk_modulus";
    info.description = pylith::topology::Field::Description(
        "bulk_modulus",
        "bulk_modulus",
        componentNames,
        componentNames.size(),
        pylith::topology::Field::SCALAR,
        _data->normalizer->getPressureScale(),
        pylith::topology::FieldQuery::validatorPositive
        );
    info.fe = pylith::topology::Field::Discretization(
        1, 2, _data->auxDim, 1, false, pylith::topology::Field::DEFAULT_BASIS, pylith::topology::Field::POLYNOMIAL_SPACE, true
        );
    info.index = 2;
    _data->subfields["bulk_modulus"] = info;

    // reference_stress
    componentNames.resize(6);
    componentNames[0] = "reference_stress_xx";
    componentNames[1] = "reference_stress_yy";
    componentNames[2] = "reference_stress_zz";
    componentNames[3] = "reference_stress_xy";
    componentNames[4] = "reference_stress_yz";
    componentNames[5] = "reference_stress_xz";
    info.description = pylith::topology::Field::Description(
        "reference_stress",
        "reference_stress",
        componentNames,
        componentNames.size(),
        pylith::topology::Field::TENSOR,
        _data->normalizer->getPressureScale()
        );
    info.fe = pylith::topology::Field::Discretization(
        2, 2, _data->auxDim, _data->auxDim, false, pylith::topology::Field::DEFAULT_BASIS, pylith::topology::Field::POLYNOMIAL_SPACE, false
        );
    info.index = 3;
    _data->subfields["reference_stress"] = info;
    if (2 == _data->auxDim) {
        _data->subfields["reference_stress"].description.numComponents = 4;
        _data->subfields["reference_stress"].description.vectorFieldType = pylith::topology::Field::OTHER;
    } // if

    // reference_strain
    componentNames.resize(6);
    componentNames[0] = "reference_strain_xx";
    componentNames[1] = "reference_strain_yy";
    componentNames[2] = "reference_strain_zz";
    componentNames[3] = "reference_strain_xy";
    componentNames[4] = "reference_strain_yz";
    componentNames[5] = "reference_strain_xz";
    info.description = pylith::topology::Field::Description(
        "reference_strain",
        "reference_strain",
        componentNames,
        componentNames.size(),
        pylith::topology::Field::TENSOR,
        1.0
        );
    info.fe = pylith::topology::Field::Discretization(
        2, 2, _data->auxDim, _data->auxDim, false, pylith::topology::Field::DEFAULT_BASIS, pylith::topology::Field::POLYNOMIAL_SPACE, false
        );
    info.index = 4;
    _data->subfields["reference_strain"] = info;
    if (2 == _data->auxDim) {
        _data->subfields["reference_strain"].description.numComponents = 4;
        _data->subfields["reference_strain"].description.vectorFieldType = pylith::topology::Field::OTHER;
    } // if

    _initialize();

    PYLITH_METHOD_END;
} // setUp


// ------------------------------------------------------------------------------------------------
// Destructor.
pylith::materials::TestAuxiliaryFactoryLinearElastic::~TestAuxiliaryFactoryLinearElastic(void) {
    PYLITH_METHOD_BEGIN;

    delete _factory;_factory = NULL;
    delete _data;_data = NULL;
    delete _mesh;_mesh = NULL;
    delete _auxiliaryField;_auxiliaryField = NULL;

    PYLITH_METHOD_END;
} // tearDown


// ------------------------------------------------------------------------------------------------
// Test adding shear modulus, bulk modulus, and reference stress/strain subfields.
void
pylith::materials::TestAuxiliaryFactoryLinearElastic::testAdd(void) {
    PYLITH_METHOD_BEGIN;
    assert(_factory);
    assert(_data);

    CHECK(!_auxiliaryField->hasSubfield("density"));
    CHECK(!_auxiliaryField->hasSubfield("shear_modulus"));
    CHECK(!_auxiliaryField->hasSubfield("bulk_modulus"));
    CHECK(!_auxiliaryField->hasSubfield("reference_stress"));
    CHECK(!_auxiliaryField->hasSubfield("reference_strain"));

    _factory->addDensity();
    _factory->addShearModulus();
    _factory->addBulkModulus();
    _factory->addReferenceStress();
    _factory->addReferenceStrain();

    assert(_data->normalizer);

    pylith::testing::FieldTester::checkSubfieldInfo(*_auxiliaryField, _data->subfields["density"]);
    pylith::testing::FieldTester::checkSubfieldInfo(*_auxiliaryField, _data->subfields["shear_modulus"]);
    pylith::testing::FieldTester::checkSubfieldInfo(*_auxiliaryField, _data->subfields["bulk_modulus"]);
    pylith::testing::FieldTester::checkSubfieldInfo(*_auxiliaryField, _data->subfields["reference_stress"]);
    pylith::testing::FieldTester::checkSubfieldInfo(*_auxiliaryField, _data->subfields["reference_strain"]);

    PYLITH_METHOD_END;
} // testAdd


// ------------------------------------------------------------------------------------------------
// Test setValues().
void
pylith::materials::TestAuxiliaryFactoryLinearElastic::testSetValuesFromDB(void) {
    PYLITH_METHOD_BEGIN;

    assert(_factory);

    _factory->addDensity();
    _factory->addShearModulus();
    _factory->addBulkModulus();
    _factory->addReferenceStress();
    _factory->addReferenceStrain();
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
pylith::materials::TestAuxiliaryFactoryLinearElastic::_initialize(void) {
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

    _factory = new AuxiliaryFactoryElastic();
    assert(_data->auxiliaryDB);
    _factory->setQueryDB(_data->auxiliaryDB);
    typedef std::map<std::string, pylith::topology::Field::SubfieldInfo>::const_iterator subfield_iter;
    for (subfield_iter iter = _data->subfields.begin(); iter != _data->subfields.end(); ++iter) {
        const char* subfieldName = iter->first.c_str();
        const pylith::topology::Field::Discretization& fe = iter->second.fe;
        _factory->setSubfieldDiscretization(subfieldName, fe.basisOrder, fe.quadOrder, fe.dimension, fe.isFaultOnly,
                                            fe.cellBasis, fe.feSpace, fe.isBasisContinuous);
    } // for
    assert(_data->normalizer);
    _factory->initialize(_auxiliaryField, *_data->normalizer, _data->dimension);

    PYLITH_METHOD_END;
} // _initialize


// ------------------------------------------------------------------------------------------------
pylith::materials::TestAuxiliaryFactoryLinearElastic_Data::TestAuxiliaryFactoryLinearElastic_Data(void) :
    meshFilename(NULL),
    cs(NULL),
    normalizer(new spatialdata::units::Nondimensional),
    auxiliaryDB(new spatialdata::spatialdb::UserFunctionDB) {}


// ------------------------------------------------------------------------------------------------
pylith::materials::TestAuxiliaryFactoryLinearElastic_Data::~TestAuxiliaryFactoryLinearElastic_Data(void) {
    delete cs;cs = NULL;
    delete normalizer;normalizer = NULL;
    delete auxiliaryDB;auxiliaryDB = NULL;
} // destructor


// End of file
