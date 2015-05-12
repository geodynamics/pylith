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

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/MeshOps.hh" // USES MeshOps::nondimensionalize()
#include "pylith/topology/Field.hh" // USES Field
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
} // setUp


// ----------------------------------------------------------------------
// Deallocate testing data.
void
pylith::materials::TestIsotropicLinearElasticityPlaneStrain::tearDown(void)
{ // tearDown
  delete _material; _material = NULL;
  delete _data; _data = NULL;
  delete _mesh; _mesh = NULL;
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
// Test preinitialize().
void
pylith::materials::TestIsotropicLinearElasticityPlaneStrain::testPreinitialize(void)
{ // testPreinitialize
  PYLITH_METHOD_BEGIN;

  CPPUNIT_ASSERT(_material);

  // Call preinitialize()
  _material->preinitialize();
  
  // Check result
  CPPUNIT_ASSERT(false); // :TODO: ADD MORE HERE

  PYLITH_METHOD_END;
} // testPreinitialize


// ----------------------------------------------------------------------
// Test _setFEKernels().
void
pylith::materials::TestIsotropicLinearElasticityPlaneStrain::test_setFEKernels(void)
{ // test_setFEKernels
  PYLITH_METHOD_BEGIN;

  CPPUNIT_ASSERT(_material);

  // Call _setFEKernels()
  // _material->_setFEKernels(solnField, prob);

  // Check result
  CPPUNIT_ASSERT(false); // :TODO: ADD MORE HERE

  PYLITH_METHOD_END;
} // test_setFEKernels


// ----------------------------------------------------------------------
// Test _dbToAuxFields().
void
pylith::materials::TestIsotropicLinearElasticityPlaneStrain::test_dbToAuxFields(void)
{ // test_dbToAuxFields
  PYLITH_METHOD_BEGIN;

  CPPUNIT_ASSERT(_material);

  // Call _dbToAuxFields()
  // :TODO: ADD MORE HERE

  // Check result
  CPPUNIT_ASSERT(false); // :TODO: ADD MORE HERE

  PYLITH_METHOD_END;
} // test_dbToAuxFields


// IntegratorPointwise ========================================


// ----------------------------------------------------------------------
// Test auxFields().
void
pylith::materials::TestIsotropicLinearElasticityPlaneStrain::testAuxFields(void)
{ // testAuxFields
  PYLITH_METHOD_BEGIN;

  CPPUNIT_ASSERT(_material);
#if 0
  // Call auxFields()
  const topology::Field& auxFields = _material->auxFields();
#endif

  // Check result
  CPPUNIT_ASSERT(false); // :TODO: ADD MORE HERE

  PYLITH_METHOD_END;
} // testAuxFields


// ----------------------------------------------------------------------
// Test hasAuxField().
void
pylith::materials::TestIsotropicLinearElasticityPlaneStrain::testHasAuxField(void)
{ // testHasAuxField
  PYLITH_METHOD_BEGIN;

  CPPUNIT_ASSERT(_material);

#if 0
  CPPUNIT_ASSERT(_material->hasAuxField("density"));
  CPPUNIT_ASSERT(_material->hasAuxField("mu"));
  CPPUNIT_ASSERT(_material->hasAuxField("lambda"));

  CPPUNIT_ASSERT(!_material->hasAuxField("abc"));
#endif

  PYLITH_METHOD_END;
} // testHaxAuxField


// ----------------------------------------------------------------------
// Test getAuxField().
void
pylith::materials::TestIsotropicLinearElasticityPlaneStrain::testGetAuxField(void)
{ // testGetAuxField
  PYLITH_METHOD_BEGIN;

  CPPUNIT_ASSERT(_material);

  // Call getAuxField()
#if 0
  topology::Field density;
  _material->getAuxField(&density, "density");
#endif

  // Check result
  CPPUNIT_ASSERT(false); // :TODO: ADD MORE HERE

  PYLITH_METHOD_END;
} // testGetAuxField


// ----------------------------------------------------------------------
// Test normalizer().
void
pylith::materials::TestIsotropicLinearElasticityPlaneStrain::testNormalizer(void)
{ // testNormalizer
  PYLITH_METHOD_BEGIN;

  CPPUNIT_ASSERT(false); // :TODO: ADD MORE HERE

  PYLITH_METHOD_END;
} // testNormalizer


// ----------------------------------------------------------------------
// Test IsJacobianSymmetric().
void
pylith::materials::TestIsotropicLinearElasticityPlaneStrain::testIsJacobianSymmetric(void)
{ // testIsJacobianSymmetric
  PYLITH_METHOD_BEGIN;

  CPPUNIT_ASSERT(_material);

  const bool flag = false; // default
  CPPUNIT_ASSERT_EQUAL(flag, _material->isJacobianSymmetric());

  // :TODO: ADD MORE HERE
  // Does flag change based on settings?

  PYLITH_METHOD_END;
} // testIsJacobianSymmetric


// ----------------------------------------------------------------------
// Test verifyConfiguration().
void
pylith::materials::TestIsotropicLinearElasticityPlaneStrain::testVerifyConfiguration(void)
{ // testVerifyConfiguration
  PYLITH_METHOD_BEGIN;

  CPPUNIT_ASSERT(_material);

  // Call verifyConfiguration()
  // :TODO: ADD MORE HERE

  // Check result
  CPPUNIT_ASSERT(false); // :TODO: ADD MORE HERE

  PYLITH_METHOD_END;
} // testVerifyConfiguration


// ----------------------------------------------------------------------
// Test initialize().
void
pylith::materials::TestIsotropicLinearElasticityPlaneStrain::testInitialize(void)
{ // testInitialize
  PYLITH_METHOD_BEGIN;

  // Call 
  CPPUNIT_ASSERT(false); // :TODO: ADD MORE HERE

  PYLITH_METHOD_END;
} // testInitialize


// ----------------------------------------------------------------------
// Test integrateResidual().
void
pylith::materials::TestIsotropicLinearElasticityPlaneStrain::testIntegrateResidual(void)
{ // testIntegrateResidual
  PYLITH_METHOD_BEGIN;

  CPPUNIT_ASSERT(false); // :TODO: ADD MORE HERE

  PYLITH_METHOD_END;
} // testIntegrateResidual


// ----------------------------------------------------------------------
// Test integrateJacobian().
void
pylith::materials::TestIsotropicLinearElasticityPlaneStrain::testIntegrateJacobian(void)
{ // testIntegrateJacobian
  PYLITH_METHOD_BEGIN;

  CPPUNIT_ASSERT(false); // :TODO: ADD MORE HERE

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
// Test dbAuxFields().
void
pylith::materials::TestIsotropicLinearElasticityPlaneStrain::testDBAuxFields(void)
{ // testDBAuxFields
  PYLITH_METHOD_BEGIN;

  CPPUNIT_ASSERT(false); // :TODO: ADD MORE HERE

  PYLITH_METHOD_END;
} // testDBAuxFields


// ----------------------------------------------------------------------
// Test _initializeAuxFieldsFromDB().
void
pylith::materials::TestIsotropicLinearElasticityPlaneStrain::test_initializeAuxFieldsFromDB(void)
{ // test_initializeAuxFieldsFromDB
  PYLITH_METHOD_BEGIN;

  CPPUNIT_ASSERT(false); // :TODO: ADD MORE HERE

  PYLITH_METHOD_END;
} // test_initializeAuxFieldsFromDB


// ----------------------------------------------------------------------
// Test _nondimAuxFields().
void
pylith::materials::TestIsotropicLinearElasticityPlaneStrain::test_nondimAuxFields(void)
{ // test_nondimAuxFields
  PYLITH_METHOD_BEGIN;

  CPPUNIT_ASSERT(false); // :TODO: ADD MORE HERE

  PYLITH_METHOD_END;
} // test_nondimAuxFields


// ----------------------------------------------------------------------
// Test _dimAuxFields().
void
pylith::materials::TestIsotropicLinearElasticityPlaneStrain::test_dimAuxFields(void)
{ // test_dimAuxFields
  PYLITH_METHOD_BEGIN;

  CPPUNIT_ASSERT(false); // :TODO: ADD MORE HERE

  PYLITH_METHOD_END;
} // test_dimAuxFields


// ----------------------------------------------------------------------
// Initialize test data.
void
pylith::materials::TestIsotropicLinearElasticityPlaneStrain::_initialize(void)
{ // _initialize
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
} // _initialize


// End of file 
