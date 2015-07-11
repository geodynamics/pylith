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

#include "TestIsotropicLinearElasticityPlaneStrain.hh" // Implementation of class methods

#include "data/IsotropicLinearElasticityPlaneStrainData.hh" // USES IsotropicLinearElasticityPlaneStrainData

#include "pylith/materials/IsotropicLinearElasticityPlaneStrain.hh" // USES IsotropicLinearElasticityPlaneStrain
#include "pylith/materials/Query.hh" // USES Query

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/MeshOps.hh" // USES MeshOps::nondimensionalize()
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/FieldQuery.hh" // USES FieldQuery
#include "pylith/meshio/MeshIOAscii.hh" // USES MeshIOAscii
#include "pylith/utils/error.hh" // USES PYLITH_METHOD_BEGIN/END

#include "spatialdata/spatialdb/SimpleDB.hh" // USES SimpleDB
#include "spatialdata/spatialdb/SimpleIOAscii.hh" // USES SimpleIOAscii
#include "spatialdata/geocoords/CSCart.hh" // USES CSCart
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional


// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::materials::TestIsotropicLinearElasticityPlaneStrain::setUp(void)
{ // setUp
  _material = new IsotropicLinearElasticityPlaneStrain();
  _data = NULL;
  _mesh = new topology::Mesh();
  _solution = NULL;
  _db = NULL;
} // setUp


// ----------------------------------------------------------------------
// Deallocate testing data.
void
pylith::materials::TestIsotropicLinearElasticityPlaneStrain::tearDown(void)
{ // tearDown
  delete _solution; _solution = NULL;
  delete _material; _material = NULL;
  delete _data; _data = NULL;
  delete _mesh; _mesh = NULL;
  delete _db; _db = NULL;
} // tearDown


// ----------------------------------------------------------------------
// Test useInertia().
void
pylith::materials::TestIsotropicLinearElasticityPlaneStrain::testUseInertia(void)
{ // testUseInertia
  PYLITH_METHOD_BEGIN;

  CPPUNIT_ASSERT(_material);

  const bool flag = false; // default
  CPPUNIT_ASSERT_EQUAL(flag, _material->_useInertia);

  _material->useInertia(!flag);
  CPPUNIT_ASSERT_EQUAL(!flag, _material->_useInertia);

  PYLITH_METHOD_END;
} // testUseInertia


// ----------------------------------------------------------------------
// Test useBodyForce().
void
pylith::materials::TestIsotropicLinearElasticityPlaneStrain::testUseBodyForce(void)
{ // testUseBodyForce
  PYLITH_METHOD_BEGIN;

  CPPUNIT_ASSERT(_material);

  const bool flag = false; // default
  CPPUNIT_ASSERT_EQUAL(flag, _material->_useBodyForce);

  _material->useBodyForce(!flag);
  CPPUNIT_ASSERT_EQUAL(!flag, _material->_useBodyForce);

  PYLITH_METHOD_END;
} // testUseBodyForce


// ----------------------------------------------------------------------
// Test auxFieldsSetup().
void
pylith::materials::TestIsotropicLinearElasticityPlaneStrain::test_auxFieldsSetup(void)
{ // test_auxFieldsSetup
  PYLITH_METHOD_BEGIN;

  CPPUNIT_ASSERT(_material);
  CPPUNIT_ASSERT(_mesh);
  CPPUNIT_ASSERT(_data);
  delete _material->_auxFields; _material->_auxFields = new topology::Field(*_mesh);CPPUNIT_ASSERT(_material->_auxFields);
  delete _material->_auxFieldsQuery; _material->_auxFieldsQuery = new topology::FieldQuery(*_material->_auxFields);CPPUNIT_ASSERT(_material->_auxFieldsQuery);
  _material->_auxFieldsSetup();
  
  // Check discretizations
  { // density
    const char* label = "density";
    const pylith::topology::Field::SubfieldInfo& info = _material->_auxFields->subfieldInfo(label);
    CPPUNIT_ASSERT_EQUAL(1, info.numComponents);
    CPPUNIT_ASSERT_EQUAL(std::string(label), info.metadata.label);
    CPPUNIT_ASSERT_EQUAL(pylith::topology::Field::SCALAR, info.metadata.vectorFieldType);
    CPPUNIT_ASSERT_EQUAL(_data->densityScale, info.metadata.scale);
    CPPUNIT_ASSERT_EQUAL(false, info.metadata.dimsOkay);
    CPPUNIT_ASSERT_EQUAL(-1, info.fe.basisOrder);
    CPPUNIT_ASSERT_EQUAL(-1, info.fe.quadOrder);
    CPPUNIT_ASSERT_EQUAL(true, info.fe.isBasisContinuous);
  } // density

  { // mu
    const char* label = "mu";
    const pylith::topology::Field::SubfieldInfo& info = _material->_auxFields->subfieldInfo(label);
    CPPUNIT_ASSERT_EQUAL(1, info.numComponents);
    CPPUNIT_ASSERT_EQUAL(std::string(label), info.metadata.label);
    CPPUNIT_ASSERT_EQUAL(pylith::topology::Field::SCALAR, info.metadata.vectorFieldType);
    CPPUNIT_ASSERT_EQUAL(_data->pressureScale, info.metadata.scale);
    CPPUNIT_ASSERT_EQUAL(false, info.metadata.dimsOkay);
    CPPUNIT_ASSERT_EQUAL(-1, info.fe.basisOrder);
    CPPUNIT_ASSERT_EQUAL(-1, info.fe.quadOrder);
    CPPUNIT_ASSERT_EQUAL(true, info.fe.isBasisContinuous);
  } // mu

  { // lambda
    const char* label = "lambda";
    const pylith::topology::Field::SubfieldInfo& info = _material->_auxFields->subfieldInfo(label);
    CPPUNIT_ASSERT_EQUAL(1, info.numComponents);
    CPPUNIT_ASSERT_EQUAL(std::string(label), info.metadata.label);
    CPPUNIT_ASSERT_EQUAL(pylith::topology::Field::SCALAR, info.metadata.vectorFieldType);
    CPPUNIT_ASSERT_EQUAL(_data->pressureScale, info.metadata.scale);
    CPPUNIT_ASSERT_EQUAL(false, info.metadata.dimsOkay);
    CPPUNIT_ASSERT_EQUAL(-1, info.fe.basisOrder);
    CPPUNIT_ASSERT_EQUAL(-1, info.fe.quadOrder);
    CPPUNIT_ASSERT_EQUAL(true, info.fe.isBasisContinuous);
  } // lambda

  if (_data->useBodyForce) { // body force
    const char* label = "body force";
    const PylithReal forceScale = _data->densityScale * _data->lengthScale / (_data->timeScale * _data->timeScale);
    const pylith::topology::Field::SubfieldInfo& info = _material->_auxFields->subfieldInfo(label);
    CPPUNIT_ASSERT_EQUAL(_data->dimension, info.numComponents);
    CPPUNIT_ASSERT_EQUAL(std::string(label), info.metadata.label);
    CPPUNIT_ASSERT_EQUAL(pylith::topology::Field::VECTOR, info.metadata.vectorFieldType);
    CPPUNIT_ASSERT_EQUAL(forceScale, info.metadata.scale);
    CPPUNIT_ASSERT_EQUAL(false, info.metadata.dimsOkay);
    CPPUNIT_ASSERT_EQUAL(-1, info.fe.basisOrder);
    CPPUNIT_ASSERT_EQUAL(-1, info.fe.quadOrder);
    CPPUNIT_ASSERT_EQUAL(true, info.fe.isBasisContinuous);
  } // body force

  // Check queries
  CPPUNIT_ASSERT_EQUAL(&Query::dbQueryDensity2D, _material->_auxFieldsQuery->queryFn("density"));
  CPPUNIT_ASSERT_EQUAL(&Query::dbQueryMu2D, _material->_auxFieldsQuery->queryFn("mu"));
  CPPUNIT_ASSERT_EQUAL(&Query::dbQueryLambda2D, _material->_auxFieldsQuery->queryFn("lambda"));
  if (_data->useBodyForce) {
    CPPUNIT_ASSERT_EQUAL(&Query::dbQueryBodyForce2D, _material->_auxFieldsQuery->queryFn("body force"));
  } // if

  PYLITH_METHOD_END;
} // test_auxFieldsSetup


// ----------------------------------------------------------------------
// Test _setFEKernels().
void
pylith::materials::TestIsotropicLinearElasticityPlaneStrain::test_setFEKernels(void)
{ // test_setFEKernels
  PYLITH_METHOD_BEGIN;

  _initializeFull(); // includes setting FE kernels

  // Check result
  CPPUNIT_ASSERT_MESSAGE("Test incomplete.", false); // :TODO: ADD MORE HERE

  PYLITH_METHOD_END;
} // test_setFEKernels


// ----------------------------------------------------------------------
// Test auxFields().
void
pylith::materials::TestIsotropicLinearElasticityPlaneStrain::testAuxFields(void)
{ // testAuxFields
  PYLITH_METHOD_BEGIN;

  CPPUNIT_ASSERT(_data);

  // Call initialize()
  _initializeFull(); // includes setting up auxFields

  CPPUNIT_ASSERT(_material);
  const pylith::topology::Field& auxFields = _material->auxFields();

  _material->_auxFields->view("AUX FIELDS"); // :TEMPORARY: Debugging

  // Check result
  CPPUNIT_ASSERT_EQUAL(std::string("auxiliary fields"), std::string(auxFields.label()));
  CPPUNIT_ASSERT_EQUAL(_data->dimension, auxFields.spaceDim());
#if 1
  PylithReal norm = 0.0;
  const PetscDM dm = auxFields.dmMesh();CPPUNIT_ASSERT(dm);
  pylith::topology::FieldQuery* query = _material->_auxFieldsQuery;
  query->openDB(_db, _data->lengthScale);
  PetscErrorCode err = DMPlexComputeL2Diff(dm, query->functions(), (void**)query->contextPtrs(), auxFields.globalVector(), &norm);PYLITH_CHECK_ERROR(err);
  query->closeDB(_db);
  const PylithReal tolerance = 1.0e-6;
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, norm, tolerance);
#else
  CPPUNIT_ASSERT_MESSAGE("Test incomplete.", false); // :TODO: ADD MORE HERE
#endif

  PYLITH_METHOD_END;
} // testAuxFields


// IntegratorPointwise ========================================


// ----------------------------------------------------------------------
// Test hasAuxField().
void
pylith::materials::TestIsotropicLinearElasticityPlaneStrain::testHasAuxField(void)
{ // testHasAuxField
  PYLITH_METHOD_BEGIN;

  _initializeFull();

  CPPUNIT_ASSERT(_material);
  CPPUNIT_ASSERT(_material->hasAuxField("density"));
  CPPUNIT_ASSERT(_material->hasAuxField("mu"));
  CPPUNIT_ASSERT(_material->hasAuxField("lambda"));

  CPPUNIT_ASSERT(!_material->hasAuxField("abc"));

  PYLITH_METHOD_END;
} // testHaxAuxField


// ----------------------------------------------------------------------
// Test getAuxField().
void
pylith::materials::TestIsotropicLinearElasticityPlaneStrain::testGetAuxField(void)
{ // testGetAuxField
  PYLITH_METHOD_BEGIN;

  _initializeFull();

  // Call getAuxField()
  CPPUNIT_ASSERT(_material);
  CPPUNIT_ASSERT(_mesh);
  pylith::topology::Field density(*_mesh);
  _material->getAuxField(&density, "density");

  density.view("DENSITY");
  // Check result
  CPPUNIT_ASSERT_MESSAGE("Test incomplete.", false); // :TODO: ADD MORE HERE

  PYLITH_METHOD_END;
} // testGetAuxField


// ----------------------------------------------------------------------
// Test normalizer().
void
pylith::materials::TestIsotropicLinearElasticityPlaneStrain::testNormalizer(void)
{ // testNormalizer
  PYLITH_METHOD_BEGIN;
  
  spatialdata::units::Nondimensional normalizer;
  const double scale = 5.0;
  normalizer.lengthScale(scale);

  CPPUNIT_ASSERT(_material);
  _material->normalizer(normalizer);
  CPPUNIT_ASSERT_EQUAL(scale, _material->_normalizer->lengthScale());

  PYLITH_METHOD_END;
} // testNormalizer


// ----------------------------------------------------------------------
// Test IsJacobianSymmetric().
void
pylith::materials::TestIsotropicLinearElasticityPlaneStrain::testIsJacobianSymmetric(void)
{ // testIsJacobianSymmetric
  PYLITH_METHOD_BEGIN;

  CPPUNIT_ASSERT(_material);

  const bool flag = true; // default
  CPPUNIT_ASSERT_EQUAL(flag, _material->isJacobianSymmetric());

  _material->useInertia(true);
  CPPUNIT_ASSERT_EQUAL(false, _material->isJacobianSymmetric());

  PYLITH_METHOD_END;
} // testIsJacobianSymmetric


// ----------------------------------------------------------------------
// Test verifyConfiguration().
void
pylith::materials::TestIsotropicLinearElasticityPlaneStrain::testVerifyConfiguration(void)
{ // testVerifyConfiguration
  PYLITH_METHOD_BEGIN;

  CPPUNIT_ASSERT(_material);

  _initializeMin();

  // Call verifyConfiguration()
  CPPUNIT_ASSERT(_material);
  CPPUNIT_ASSERT(_mesh);
  _material->verifyConfiguration(*_mesh);

  // Nothing to test.

  PYLITH_METHOD_END;
} // testVerifyConfiguration


// ----------------------------------------------------------------------
// Test integrateResidual().
void
pylith::materials::TestIsotropicLinearElasticityPlaneStrain::testIntegrateResidual(void)
{ // testIntegrateResidual
  PYLITH_METHOD_BEGIN;

  CPPUNIT_ASSERT_MESSAGE("Test incomplete.", false); // :TODO: ADD MORE HERE

  PYLITH_METHOD_END;
} // testIntegrateResidual


// ----------------------------------------------------------------------
// Test integrateJacobian().
void
pylith::materials::TestIsotropicLinearElasticityPlaneStrain::testIntegrateJacobian(void)
{ // testIntegrateJacobian
  PYLITH_METHOD_BEGIN;

  CPPUNIT_ASSERT_MESSAGE("Test incomplete.", false); // :TODO: ADD MORE HERE

  PYLITH_METHOD_END;
} // testIntegrateJacobian

// MaterialNew ========================================

// ----------------------------------------------------------------------
// Test dimension().
void
pylith::materials::TestIsotropicLinearElasticityPlaneStrain::testDimension(void)
{ // testDimension
  PYLITH_METHOD_BEGIN;

  CPPUNIT_ASSERT(_material);
  CPPUNIT_ASSERT(_data);
  CPPUNIT_ASSERT_EQUAL(_data->dimension, _material->dimension());

  PYLITH_METHOD_END;
} // testDimension


// ----------------------------------------------------------------------
// Test id().
void
pylith::materials::TestIsotropicLinearElasticityPlaneStrain::testId(void)
{ // testId
  PYLITH_METHOD_BEGIN;

  CPPUNIT_ASSERT(_material);
  CPPUNIT_ASSERT(_data);
  CPPUNIT_ASSERT_EQUAL(_data->id, _material->id());

  PYLITH_METHOD_END;
} // testId


// ----------------------------------------------------------------------
// Test label().
void
pylith::materials::TestIsotropicLinearElasticityPlaneStrain::testLabel(void)
{ // testLabel
  PYLITH_METHOD_BEGIN;

  CPPUNIT_ASSERT(_material);
  CPPUNIT_ASSERT(_data);
  CPPUNIT_ASSERT_EQUAL(std::string(_data->label), std::string(_material->label()));

  PYLITH_METHOD_END;
} // testLabel


// ----------------------------------------------------------------------
// Test auxFieldsDB().
void
pylith::materials::TestIsotropicLinearElasticityPlaneStrain::testAuxFieldsDB(void)
{ // testAuxFieldsDB
  PYLITH_METHOD_BEGIN;

  const std::string label = "test db";
  spatialdata::spatialdb::SimpleDB db;
  db.label(label.c_str());

  CPPUNIT_ASSERT(_material);
  _material->auxFieldsDB(&db);

  CPPUNIT_ASSERT(_material->_auxFieldsDB);
  CPPUNIT_ASSERT_EQUAL(label, std::string(_material->_auxFieldsDB->label()));

  PYLITH_METHOD_END;
} // testAuxFieldsDB


// ----------------------------------------------------------------------
// Test discretization().
void
pylith::materials::TestIsotropicLinearElasticityPlaneStrain::testDiscretization(void)
{ // testDiscretization
  PYLITH_METHOD_BEGIN;

  const topology::FieldBase::DiscretizeInfo infoDefault = {-1, -1, true};
  const topology::FieldBase::DiscretizeInfo infoA = {1, 2, false};
  const topology::FieldBase::DiscretizeInfo infoB = {2, 2, true};
  
  CPPUNIT_ASSERT(_material);
  _material->discretization("A", infoA);
  _material->discretization("B", infoB);

  { // A
    const topology::FieldBase::DiscretizeInfo& test = _material->discretization("A");
    CPPUNIT_ASSERT_EQUAL(test.basisOrder, infoA.basisOrder);
    CPPUNIT_ASSERT_EQUAL(test.quadOrder, infoA.quadOrder);
    CPPUNIT_ASSERT_EQUAL(test.isBasisContinuous, infoA.isBasisContinuous);
  } // A

  { // B
    const topology::FieldBase::DiscretizeInfo& test = _material->discretization("B");
    CPPUNIT_ASSERT_EQUAL(test.basisOrder, infoB.basisOrder);
    CPPUNIT_ASSERT_EQUAL(test.quadOrder, infoB.quadOrder);
    CPPUNIT_ASSERT_EQUAL(test.isBasisContinuous, infoB.isBasisContinuous);
  } // B

  { // C (default)
    const topology::FieldBase::DiscretizeInfo& test = _material->discretization("C");
    CPPUNIT_ASSERT_EQUAL(test.basisOrder, infoDefault.basisOrder);
    CPPUNIT_ASSERT_EQUAL(test.quadOrder, infoDefault.quadOrder);
    CPPUNIT_ASSERT_EQUAL(test.isBasisContinuous, infoDefault.isBasisContinuous);
  } // C (default)

  { // default
    const topology::FieldBase::DiscretizeInfo& test = _material->discretization("default");
    CPPUNIT_ASSERT_EQUAL(test.basisOrder, infoDefault.basisOrder);
    CPPUNIT_ASSERT_EQUAL(test.quadOrder, infoDefault.quadOrder);
    CPPUNIT_ASSERT_EQUAL(test.isBasisContinuous, infoDefault.isBasisContinuous);
  } // default

  PYLITH_METHOD_END;
} // testDiscretization


// ----------------------------------------------------------------------
// Initialize test data.
void
pylith::materials::TestIsotropicLinearElasticityPlaneStrain::_initializeMin(void)
{ // _initializeMin
  PYLITH_METHOD_BEGIN;

  CPPUNIT_ASSERT(_material);
  CPPUNIT_ASSERT(_data);

  meshio::MeshIOAscii iohandler;
  iohandler.filename(_data->filenameMesh);
  iohandler.read(_mesh);CPPUNIT_ASSERT(_mesh);

  // Setup coordinates.
  spatialdata::geocoords::CSCart cs;
  cs.setSpaceDim(_mesh->dimension());
  cs.initialize();
  _mesh->coordsys(&cs);

  // Setup scales.
  spatialdata::units::Nondimensional normalizer;
  normalizer.lengthScale(_data->lengthScale);
  normalizer.pressureScale(_data->pressureScale);
  normalizer.timeScale(_data->timeScale);
  normalizer.densityScale(_data->densityScale);
  topology::MeshOps::nondimensionalize(_mesh, normalizer);

  _material->id(_data->id);
  _material->label(_data->label);
  _material->useInertia(_data->useInertia);
  _material->useBodyForce(_data->useBodyForce);
  _material->normalizer(normalizer);

  PYLITH_METHOD_END;
} // _initializeMin


// ----------------------------------------------------------------------
// Initialize test data.
void
pylith::materials::TestIsotropicLinearElasticityPlaneStrain::_initializeFull(void)
{ // _initializeFull
  PYLITH_METHOD_BEGIN;

  CPPUNIT_ASSERT(_material);
  CPPUNIT_ASSERT(_data);
  CPPUNIT_ASSERT(_mesh);

  _initializeMin();

  // Set auxiliary fields spatial database.
  delete _db; _db = new spatialdata::spatialdb::SimpleDB;CPPUNIT_ASSERT(_db);
  spatialdata::spatialdb::SimpleIOAscii dbIO;
  dbIO.filename(_data->filenameAuxFieldsDB);
  _db->ioHandler(&dbIO);
  _db->label("IsotropicLinearElasciticityPlaneStrain auxiliary fields");
  _db->queryType(spatialdata::spatialdb::SimpleDB::NEAREST);
  _material->auxFieldsDB(_db);

  // Create solution field.
  delete _solution; _solution = new pylith::topology::Field(*_mesh);
  CPPUNIT_ASSERT(_data->discretizations);
  _solution->subfieldAdd("displacement", _data->dimension, topology::Field::VECTOR, _data->discretizations[0], _data->lengthScale);
  if (_data->useInertia) {
    const PylithReal velocityScale = _data->lengthScale / _data->timeScale;
    _solution->subfieldAdd("velocity", _data->dimension, topology::Field::VECTOR, _data->discretizations[1], velocityScale);
  } // if
  _solution->subfieldsSetup();
  _solution->allocate();
  _solution->zeroAll();

  CPPUNIT_ASSERT(_solution);
  _material->initialize(*_solution);

  PYLITH_METHOD_END;
} // _initializeFull


// End of file 
