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
#include "pylith/topology/VisitorMesh.hh" // USES VecVisitorMesh
#include "pylith/topology/FieldQuery.hh" // USES FieldQuery
#include "pylith/topology/Jacobian.hh" // USES Jacobian
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
  _solution1 = NULL;
  _solution2 = NULL;
  _solution1Dot = NULL;
  _solution2Dot = NULL;
  _auxDB = NULL;
  _soln1DB = NULL;
  _soln2DB = NULL;
} // setUp


// ----------------------------------------------------------------------
// Deallocate testing data.
void
pylith::materials::TestIsotropicLinearElasticityPlaneStrain::tearDown(void)
{ // tearDown
  delete _solution1; _solution1 = NULL;
  delete _solution2; _solution2 = NULL;
  delete _solution1Dot; _solution1Dot = NULL;
  delete _solution2Dot; _solution2Dot = NULL;
  delete _material; _material = NULL;
  delete _data; _data = NULL;
  delete _mesh; _mesh = NULL;
  delete _auxDB; _auxDB = NULL;
  delete _soln1DB; _soln1DB = NULL;
  delete _soln2DB; _soln2DB = NULL;
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

  { // shearModulus
    const char* label = "shaer_modulus";
    const pylith::topology::Field::SubfieldInfo& info = _material->_auxFields->subfieldInfo(label);
    CPPUNIT_ASSERT_EQUAL(1, info.numComponents);
    CPPUNIT_ASSERT_EQUAL(std::string(label), info.metadata.label);
    CPPUNIT_ASSERT_EQUAL(pylith::topology::Field::SCALAR, info.metadata.vectorFieldType);
    CPPUNIT_ASSERT_EQUAL(_data->pressureScale, info.metadata.scale);
    CPPUNIT_ASSERT_EQUAL(false, info.metadata.dimsOkay);
    CPPUNIT_ASSERT_EQUAL(-1, info.fe.basisOrder);
    CPPUNIT_ASSERT_EQUAL(-1, info.fe.quadOrder);
    CPPUNIT_ASSERT_EQUAL(true, info.fe.isBasisContinuous);
  } // shearModulus

  { // bulkModulus
    const char* label = "bulk_modulus";
    const pylith::topology::Field::SubfieldInfo& info = _material->_auxFields->subfieldInfo(label);
    CPPUNIT_ASSERT_EQUAL(1, info.numComponents);
    CPPUNIT_ASSERT_EQUAL(std::string(label), info.metadata.label);
    CPPUNIT_ASSERT_EQUAL(pylith::topology::Field::SCALAR, info.metadata.vectorFieldType);
    CPPUNIT_ASSERT_EQUAL(_data->pressureScale, info.metadata.scale);
    CPPUNIT_ASSERT_EQUAL(false, info.metadata.dimsOkay);
    CPPUNIT_ASSERT_EQUAL(-1, info.fe.basisOrder);
    CPPUNIT_ASSERT_EQUAL(-1, info.fe.quadOrder);
    CPPUNIT_ASSERT_EQUAL(true, info.fe.isBasisContinuous);
  } // bulkModulus

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

  // Make sure DB query functions are set correctly.
  CPPUNIT_ASSERT_EQUAL(&pylith::topology::FieldQuery::dbQueryGeneric, _material->_auxFieldsQuery->queryFn("density"));
  CPPUNIT_ASSERT_EQUAL(&pylith::materials::Query::dbQueryShearModulus2D, _material->_auxFieldsQuery->queryFn("shear_modulus"));
  CPPUNIT_ASSERT_EQUAL(&pylith::materials::Query::dbQueryBulkModulus2D, _material->_auxFieldsQuery->queryFn("bulk_modulus"));
  if (_data->useBodyForce) {
    CPPUNIT_ASSERT_EQUAL(&pylith::topology::FieldQuery::dbQueryGeneric, _material->_auxFieldsQuery->queryFn("body force"));
  } // if

  PYLITH_METHOD_END;
} // test_auxFieldsSetup


// ----------------------------------------------------------------------
// Test _setFEKernelsRHSResidual().
void
pylith::materials::TestIsotropicLinearElasticityPlaneStrain::test_setFEKernelsRHSResidual(void)
{ // test_setFEKernelsRHSResidual
  PYLITH_METHOD_BEGIN;

  _initializeFull();

  CPPUNIT_ASSERT(_solution1);
  CPPUNIT_ASSERT(_material);
  _material->_setFEKernelsRHSResidual(*_solution1);

  PetscDS prob = NULL;
  PetscErrorCode err = DMGetDS(_solution1->dmMesh(), &prob);CPPUNIT_ASSERT(!err);CPPUNIT_ASSERT(prob);

  const int numSolnFields = _data->numSolnFields;
  const int numFieldKernels = _data->numKernelsResidual;
  for (int iField=0; iField < numSolnFields; ++iField) {
    const int indexK = iField*numFieldKernels;

    const PetscPointFunc* kernelsE = _data->kernelsRHSResidual;CPPUNIT_ASSERT(kernelsE);
    PetscPointFunc g0 = NULL;
    PetscPointFunc g1 = NULL;
    err = PetscDSGetResidual(prob, iField, &g0, &g1);CPPUNIT_ASSERT(!err);
    CPPUNIT_ASSERT(kernelsE[indexK+0] == g0);
    CPPUNIT_ASSERT(kernelsE[indexK+1] == g1);
  } // for

  PYLITH_METHOD_END;
} // test_setFEKernelsRHSResidual


// ----------------------------------------------------------------------
// Test _setFEKernelsRHSJacobian().
void
pylith::materials::TestIsotropicLinearElasticityPlaneStrain::test_setFEKernelsRHSJacobian(void)
{ // test_setFEKernelsRHSJacobian
  PYLITH_METHOD_BEGIN;

  _initializeFull();

  CPPUNIT_ASSERT(_solution1);
  CPPUNIT_ASSERT(_material);
  _material->_setFEKernelsRHSJacobian(*_solution1);

  PetscDS prob = NULL;
  PetscErrorCode err = DMGetDS(_solution1->dmMesh(), &prob);CPPUNIT_ASSERT(!err);CPPUNIT_ASSERT(prob);

  const int numSolnFields = _data->numSolnFields;
  const int numFieldKernels = _data->numKernelsJacobian;
  for (int iField=0; iField < numSolnFields; ++iField) {
    for (int jField=0; jField < numSolnFields; ++jField) {
      const int indexK = iField*numSolnFields*numFieldKernels + jField*numFieldKernels;
 
      const PetscPointJac* kernelsE = _data->kernelsRHSJacobian;CPPUNIT_ASSERT(kernelsE);
      PetscPointJac Jg0 = NULL;
      PetscPointJac Jg1 = NULL;
      PetscPointJac Jg2 = NULL;
      PetscPointJac Jg3 = NULL;
      err = PetscDSGetJacobian(prob, iField, jField, &Jg0, &Jg1, &Jg2, &Jg3);CPPUNIT_ASSERT(!err);
      CPPUNIT_ASSERT(kernelsE[indexK+0] == Jg0);
      CPPUNIT_ASSERT(kernelsE[indexK+1] == Jg1);
      CPPUNIT_ASSERT(kernelsE[indexK+2] == Jg2);
      CPPUNIT_ASSERT(kernelsE[indexK+3] == Jg3);
    } // for
  } // for
  
  PYLITH_METHOD_END;
} // test_setFEKernelsRHSJacobian


// ----------------------------------------------------------------------
// Test _setFEKernelsLHSResidual().
void
pylith::materials::TestIsotropicLinearElasticityPlaneStrain::test_setFEKernelsLHSResidual(void)
{ // test_setFEKernelsLHSResidual
  PYLITH_METHOD_BEGIN;

  _initializeFull();

  CPPUNIT_ASSERT(_solution1);
  CPPUNIT_ASSERT(_material);
  _material->_setFEKernelsLHSResidual(*_solution1);

  PetscDS prob = NULL;
  PetscErrorCode err = DMGetDS(_solution1->dmMesh(), &prob);CPPUNIT_ASSERT(!err);CPPUNIT_ASSERT(prob);

  const int numSolnFields = _data->numSolnFields;
  const int numFieldKernels = _data->numKernelsResidual;
  for (int iField=0; iField < numSolnFields; ++iField) {
    const int indexK = iField*numFieldKernels;

    const PetscPointFunc* kernelsE = _data->kernelsLHSResidual;CPPUNIT_ASSERT(kernelsE);
    PetscPointFunc f0 = NULL;
    PetscPointFunc f1 = NULL;
    err = PetscDSGetResidual(prob, iField, &f0, &f1);CPPUNIT_ASSERT(!err);
    CPPUNIT_ASSERT(kernelsE[indexK+0] == f0);
    CPPUNIT_ASSERT(kernelsE[indexK+1] == f1);
  } // for

  PYLITH_METHOD_END;
} // test_setFEKernelsLHSResidual


// ----------------------------------------------------------------------
// Test _setFEKernelsLHSJacobianImplicit().
void
pylith::materials::TestIsotropicLinearElasticityPlaneStrain::test_setFEKernelsLHSJacobianImplicit(void)
{ // test_setFEKernelsLHSJacobianImplicit
  PYLITH_METHOD_BEGIN;

  _initializeFull();

  CPPUNIT_ASSERT(_solution1);
  CPPUNIT_ASSERT(_material);
  _material->_setFEKernelsLHSJacobianImplicit(*_solution1);

  PetscDS prob = NULL;
  PetscErrorCode err = DMGetDS(_solution1->dmMesh(), &prob);CPPUNIT_ASSERT(!err);CPPUNIT_ASSERT(prob);

  const int numSolnFields = _data->numSolnFields;
  const int numFieldKernels = _data->numKernelsJacobian;
  for (int iField=0; iField < numSolnFields; ++iField) {
    for (int jField=0; jField < numSolnFields; ++jField) {
      const int indexK = iField*numSolnFields*numFieldKernels + jField*numFieldKernels;

      const PetscPointJac* kernelsE = _data->kernelsLHSJacobianImplicit;CPPUNIT_ASSERT(kernelsE);
      PetscPointJac Jf0 = NULL;
      PetscPointJac Jf1 = NULL;
      PetscPointJac Jf2 = NULL;
      PetscPointJac Jf3 = NULL;
      err = PetscDSGetJacobian(prob, iField, jField, &Jf0, &Jf1, &Jf2, &Jf3);CPPUNIT_ASSERT(!err);
      CPPUNIT_ASSERT(kernelsE[indexK+0] == Jf0);
      CPPUNIT_ASSERT(kernelsE[indexK+1] == Jf1);
      CPPUNIT_ASSERT(kernelsE[indexK+2] == Jf2);
      CPPUNIT_ASSERT(kernelsE[indexK+3] == Jf3);
    } // for
  } // for
  
  PYLITH_METHOD_END;
} // test_setFEKernelsLHSJacobianImplicit


// ----------------------------------------------------------------------
// Test _setFEKernelsLHSJacobianExplicit().
void
pylith::materials::TestIsotropicLinearElasticityPlaneStrain::test_setFEKernelsLHSJacobianExplicit(void)
{ // test_setFEKernelsLHSJacobianExplicit
  PYLITH_METHOD_BEGIN;

  _initializeFull();

  CPPUNIT_ASSERT(_solution1);
  CPPUNIT_ASSERT(_material);
  _material->_setFEKernelsLHSJacobianExplicit(*_solution1);

  PetscDS prob = NULL;
  PetscErrorCode err = DMGetDS(_solution1->dmMesh(), &prob);CPPUNIT_ASSERT(!err);CPPUNIT_ASSERT(prob);

  const int numSolnFields = _data->numSolnFields;
  const int numFieldKernels = _data->numKernelsJacobian;
  for (int iField=0; iField < numSolnFields; ++iField) {
    for (int jField=0; jField < numSolnFields; ++jField) {
      const int indexK = iField*numSolnFields*numFieldKernels + jField*numFieldKernels;

      const PetscPointJac* kernelsE = _data->kernelsLHSJacobianExplicit;CPPUNIT_ASSERT(kernelsE);
      PetscPointJac Jf0 = NULL;
      PetscPointJac Jf1 = NULL;
      PetscPointJac Jf2 = NULL;
      PetscPointJac Jf3 = NULL;
      err = PetscDSGetJacobian(prob, iField, jField, &Jf0, &Jf1, &Jf2, &Jf3);CPPUNIT_ASSERT(!err);
      CPPUNIT_ASSERT(kernelsE[indexK+0] == Jf0);
      CPPUNIT_ASSERT(kernelsE[indexK+1] == Jf1);
      CPPUNIT_ASSERT(kernelsE[indexK+2] == Jf2);
      CPPUNIT_ASSERT(kernelsE[indexK+3] == Jf3);
    } // for
  } // for
  
  PYLITH_METHOD_END;
} // test_setFEKernelsLHSJacobianExplicit


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
  CPPUNIT_ASSERT(_material->hasAuxField("shear_modulus"));
  CPPUNIT_ASSERT(_material->hasAuxField("bulk_modulus"));

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

  CPPUNIT_ASSERT(_material);
  CPPUNIT_ASSERT(_mesh);

  { // Test getting density field.
    pylith::topology::Field density(*_mesh);
    _material->getAuxField(&density, "density");

    density.createScatter(density.mesh()); // Populate global vector.
    density.scatterLocalToGlobal();

    //density.view("DENSITY"); // DEBUGGING

    // Check result
    CPPUNIT_ASSERT_EQUAL(std::string("density"), std::string(density.label()));
    CPPUNIT_ASSERT_EQUAL(_data->dimension, density.spaceDim());

    pylith::topology::FieldQuery queryDensity(density);
    queryDensity.queryFn("density", pylith::topology::FieldQuery::dbQueryGeneric);
    queryDensity.openDB(_auxDB, _data->lengthScale);

    PylithReal norm = 0.0;
    const PylithReal t = _data->t1;
    const PetscDM dm = density.dmMesh();CPPUNIT_ASSERT(dm);
    PetscErrorCode err = DMComputeL2Diff(dm, t, queryDensity.functions(), (void**)queryDensity.contextPtrs(), density.globalVector(), &norm);CPPUNIT_ASSERT(!err);
    queryDensity.closeDB(_auxDB);

    const PylithReal tolerance = 1.0e-6;
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, norm, tolerance);
  } // Test getting density field

  { // Test getting bulkModulus field.
    pylith::topology::Field bulkModulus(*_mesh);
    _material->getAuxField(&bulkModulus, "bulk_modulus");

    bulkModulus.createScatter(bulkModulus.mesh()); // Populate global vector.
    bulkModulus.scatterLocalToGlobal();

    //bulkModulus.view("BULK MODULUS"); // DEBUGGING

    // Check result
    CPPUNIT_ASSERT_EQUAL(std::string("bulk_modulus"), std::string(bulkModulus.label()));
    CPPUNIT_ASSERT_EQUAL(_data->dimension, bulkModulus.spaceDim());

    pylith::topology::FieldQuery queryBulkModulus(bulkModulus);
    queryBulkModulus.queryFn("bulk_modulus", pylith::materials::Query::dbQueryBulkModulus2D);
    queryBulkModulus.openDB(_auxDB, _data->lengthScale);

    PylithReal norm = 0.0;
    const PylithReal t = _data->t1;
    const PetscDM dm = bulkModulus.dmMesh();CPPUNIT_ASSERT(dm);
    PetscErrorCode err = DMComputeL2Diff(dm, t, queryBulkModulus.functions(), (void**)queryBulkModulus.contextPtrs(), bulkModulus.globalVector(), &norm);CPPUNIT_ASSERT(!err);
    queryBulkModulus.closeDB(_auxDB);

    const PylithReal tolerance = 1.0e-6;
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, norm, tolerance);
  } // Test getting bulkModulus field


  PYLITH_METHOD_END;
} // testGetAuxField


// ----------------------------------------------------------------------
// Test auxFieldDiscretization().
void
pylith::materials::TestIsotropicLinearElasticityPlaneStrain::testAuxFieldsDiscretization(void)
{ // testAuxFieldsDiscretization
  PYLITH_METHOD_BEGIN;

  const topology::FieldBase::DiscretizeInfo infoDefault = {-1, -1, true};
  const topology::FieldBase::DiscretizeInfo infoA = {1, 2, false};
  const topology::FieldBase::DiscretizeInfo infoB = {2, 2, true};
  
  CPPUNIT_ASSERT(_material);
  _material->auxFieldDiscretization("A", infoA);
  _material->auxFieldDiscretization("B", infoB);

  { // A
    const topology::FieldBase::DiscretizeInfo& test = _material->auxFieldDiscretization("A");
    CPPUNIT_ASSERT_EQUAL(test.basisOrder, infoA.basisOrder);
    CPPUNIT_ASSERT_EQUAL(test.quadOrder, infoA.quadOrder);
    CPPUNIT_ASSERT_EQUAL(test.isBasisContinuous, infoA.isBasisContinuous);
  } // A

  { // B
    const topology::FieldBase::DiscretizeInfo& test = _material->auxFieldDiscretization("B");
    CPPUNIT_ASSERT_EQUAL(test.basisOrder, infoB.basisOrder);
    CPPUNIT_ASSERT_EQUAL(test.quadOrder, infoB.quadOrder);
    CPPUNIT_ASSERT_EQUAL(test.isBasisContinuous, infoB.isBasisContinuous);
  } // B

  { // C (default)
    const topology::FieldBase::DiscretizeInfo& test = _material->auxFieldDiscretization("C");
    CPPUNIT_ASSERT_EQUAL(test.basisOrder, infoDefault.basisOrder);
    CPPUNIT_ASSERT_EQUAL(test.quadOrder, infoDefault.quadOrder);
    CPPUNIT_ASSERT_EQUAL(test.isBasisContinuous, infoDefault.isBasisContinuous);
  } // C (default)

  { // default
    const topology::FieldBase::DiscretizeInfo& test = _material->auxFieldDiscretization("default");
    CPPUNIT_ASSERT_EQUAL(test.basisOrder, infoDefault.basisOrder);
    CPPUNIT_ASSERT_EQUAL(test.quadOrder, infoDefault.quadOrder);
    CPPUNIT_ASSERT_EQUAL(test.isBasisContinuous, infoDefault.isBasisContinuous);
  } // default

  PYLITH_METHOD_END;
} // testAuxFieldsDiscretization


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
// Test checkConstraints().
void
pylith::materials::TestIsotropicLinearElasticityPlaneStrain::testCheckConstraints(void)
{ // testCheckConstraints
  PYLITH_METHOD_BEGIN;

  CPPUNIT_ASSERT(_material);

  _initializeFull();

  // Call checkConstraints()
  CPPUNIT_ASSERT(_material);
  CPPUNIT_ASSERT(_mesh);
  _material->checkConstraints(*_mesh);

  // Nothing to test.

  PYLITH_METHOD_END;
} // testVerifyConfiguration


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
  CPPUNIT_ASSERT_EQUAL(_data->materialId, _material->id());

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
  CPPUNIT_ASSERT_EQUAL(std::string(_data->materialLabel), std::string(_material->label()));

  PYLITH_METHOD_END;
} // testLabel


// ----------------------------------------------------------------------
// Test initialize().
void
pylith::materials::TestIsotropicLinearElasticityPlaneStrain::testInitialize(void)
{ // testInitialize
  PYLITH_METHOD_BEGIN;

  // Call initialize()
  _initializeFull(); // includes setting up auxFields

  CPPUNIT_ASSERT(_material);
  const pylith::topology::Field& auxFields = _material->auxFields();

  //_material->_auxFields->view("AUX FIELDS"); // :DEBUGGING:

  // Check result
  CPPUNIT_ASSERT(_data);
  CPPUNIT_ASSERT_EQUAL(std::string("auxiliary fields"), std::string(auxFields.label()));
  CPPUNIT_ASSERT_EQUAL(_data->dimension, auxFields.spaceDim());

  PylithReal norm = 0.0;
  PylithReal t = 0.0;
  const PetscDM dm = auxFields.dmMesh();CPPUNIT_ASSERT(dm);
  pylith::topology::FieldQuery* query = _material->_auxFieldsQuery;
  query->openDB(_auxDB, _data->lengthScale);

  PetscErrorCode err = DMComputeL2Diff(dm, t, query->functions(), (void**)query->contextPtrs(), auxFields.globalVector(), &norm);CPPUNIT_ASSERT(!err);
  query->closeDB(_auxDB);
  const PylithReal tolerance = 1.0e-6;
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, norm, tolerance);

  PYLITH_METHOD_END;
} // testInitialize


// ----------------------------------------------------------------------
// Test computeRHSResidual().
void
pylith::materials::TestIsotropicLinearElasticityPlaneStrain::testComputeResidual(void)
{ // testComputeResidual
  PYLITH_METHOD_BEGIN;

  // Call initialize()
  _initializeFull(); // includes setting up auxFields

  CPPUNIT_ASSERT(_mesh);
  pylith::topology::Field residualRHS(*_mesh);
  residualRHS.cloneSection(*_solution1);
  residualRHS.label("residual RHS");
  residualRHS.allocate();
  residualRHS.zeroAll();

  pylith::topology::Field residualLHS(*_mesh);
  residualLHS.cloneSection(*_solution1);
  residualLHS.label("residual LHS");
  residualLHS.allocate();
  residualLHS.zeroAll();

  CPPUNIT_ASSERT(_material);
  const PylithReal t = _data->t1;
  const PylithReal dt = _data->dt;
  _material->computeRHSResidual(&residualRHS, t, dt, *_solution1);
  _material->computeLHSResidual(&residualLHS, t, dt, *_solution1, *_solution1Dot);

  _zeroBoundary(&residualRHS);
  _zeroBoundary(&residualLHS);

  // Scatter local to global.
  residualRHS.complete();
  residualLHS.complete();

  PetscVec residualVec = NULL;
  PetscErrorCode err;
  err = VecDuplicate(residualRHS.globalVector(), &residualVec);CPPUNIT_ASSERT(!err);
  err = VecZeroEntries(residualVec);CPPUNIT_ASSERT(!err);
  err = VecWAXPY(residualVec, -1.0, residualLHS.globalVector(), residualRHS.globalVector());CPPUNIT_ASSERT(!err);
  
  //residualRHS.view("RESIDUAL RHS"); // DEBUGGING
  //residualLHS.view("RESIDUAL LHS"); // DEBUGGING

  PylithReal norm = 0.0;
  err = VecNorm(residualVec, NORM_2, &norm);CPPUNIT_ASSERT(!err);
  err = VecDestroy(&residualVec);CPPUNIT_ASSERT(!err);
  const PylithReal tolerance = 1.0e-6;
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, norm, tolerance);
  CPPUNIT_ASSERT(norm > 0.0); // Norm exactly equal to zero almost certainly means test is satisfied trivially.
  
  PYLITH_METHOD_END;
} // testComputeResidual


// ----------------------------------------------------------------------
// Test computeRHSJacobian().
void
pylith::materials::TestIsotropicLinearElasticityPlaneStrain::testComputeRHSJacobian(void)
{ // testComputeRHSJacobian
  PYLITH_METHOD_BEGIN;

  /*
    Create linear problem (MMS) with two solutions, s_1 and s_2.
    Check that Jg(s_1)*(s_2 - s_1) = G(s_2) - G(s_1).
  */
  
  // Call initialize()
  _initializeFull(); // includes setting up auxFields

  CPPUNIT_ASSERT(_mesh);
  pylith::topology::Field residual1(*_mesh);
  residual1.cloneSection(*_solution1);
  residual1.label("residual1");
  residual1.allocate();
  residual1.zeroAll();

  pylith::topology::Field residual2(*_mesh);
  residual2.cloneSection(*_solution2);
  residual2.label("residual2");
  residual2.allocate();
  residual2.zeroAll();

#if 0 // DEBUGGING
  PetscOptionsSetValue(NULL, "-dm_plex_print_fem", "2");
  DMSetFromOptions(_solution1->dmMesh());
#endif

  CPPUNIT_ASSERT(_material);
  const PylithReal t1 = _data->t1;
  const PylithReal t2 = _data->t2;
  const PylithReal dt = _data->dt;
  _material->computeRHSResidual(&residual1, t1, dt, *_solution1);
  _material->computeRHSResidual(&residual2, t2, dt, *_solution2);

  // Scatter local to global.
  _solution1->createScatter(_solution1->mesh());
  _solution2->createScatter(_solution2->mesh());
  _solution1->scatterLocalToGlobal();
  _solution2->scatterLocalToGlobal();
  residual1.complete();
  residual2.complete();

  // Check that J(s_1)*(s_2 - s_1) = G(s_2) - G(s_1).

  PetscVec residualVec = NULL;
  PetscErrorCode err;
  err = VecDuplicate(residual1.globalVector(), &residualVec);CPPUNIT_ASSERT(!err);
  err = VecZeroEntries(residualVec);CPPUNIT_ASSERT(!err);
  err = VecWAXPY(residualVec, -1.0, residual1.globalVector(), residual2.globalVector());CPPUNIT_ASSERT(!err);

  PetscVec solnIncrVec = NULL;
  err = VecDuplicate(_solution1->globalVector(), &solnIncrVec);CPPUNIT_ASSERT(!err);
  err = VecZeroEntries(solnIncrVec);CPPUNIT_ASSERT(!err);
  err = VecWAXPY(solnIncrVec, -1.0, _solution1->globalVector(), _solution2->globalVector());CPPUNIT_ASSERT(!err);
  
  //residual1.view("RESIDUAL RHS"); // DEBUGGING
  //residual2.view("RESIDUAL LHS"); // DEBUGGING

  // Compute Jacobian
  pylith::topology::Jacobian jacobian(*_solution1);
  pylith::topology::Jacobian* preconditioner = &jacobian; // Use Jacobian == preconditioner.
  _material->computeRHSJacobian(&jacobian, preconditioner, t1, dt, *_solution1);
  CPPUNIT_ASSERT_EQUAL(false, _material->needNewJacobian());
  jacobian.assemble("final_assembly");

  // result = J*(-solnIncr) + residual
  PetscVec resultVec = NULL;
  err = VecDuplicate(_solution1->globalVector(), &resultVec);CPPUNIT_ASSERT(!err);
  err = VecZeroEntries(resultVec);CPPUNIT_ASSERT(!err);
  err = VecScale(solnIncrVec, -1.0);CPPUNIT_ASSERT(!err);
  err = MatMultAdd(jacobian.matrix(), solnIncrVec, residualVec, resultVec);CPPUNIT_ASSERT(!err);

#if 0 // DEBUGGING  
  std::cout << "SOLN INCR" << std::endl;
  VecView(solnIncrVec, PETSC_VIEWER_STDOUT_SELF);
  std::cout << "F2-F1" << std::endl;
  VecView(residualVec, PETSC_VIEWER_STDOUT_SELF);
  std::cout << "RESULT" << std::endl;
  VecView(resultVec, PETSC_VIEWER_STDOUT_SELF);
#endif

  PylithReal norm = 0.0;
  err = VecNorm(resultVec, NORM_2, &norm);CPPUNIT_ASSERT(!err);
  err = VecDestroy(&resultVec);CPPUNIT_ASSERT(!err);
  err = VecDestroy(&solnIncrVec);CPPUNIT_ASSERT(!err);
  err = VecDestroy(&residualVec);CPPUNIT_ASSERT(!err);

  const PylithReal tolerance = 1.0e-6;
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, norm, tolerance);
  CPPUNIT_ASSERT(norm > 0.0); // Norm exactly equal to zero almost certainly means test is satisfied trivially.
  
  PYLITH_METHOD_END;
} // testComputeRHSJacobian


// ----------------------------------------------------------------------
// Test computeLHSJacobianImplicit().
void
pylith::materials::TestIsotropicLinearElasticityPlaneStrain::testComputeLHSJacobianImplicit(void)
{ // testComputeLHSJacobianImplicit
  PYLITH_METHOD_BEGIN;

  CPPUNIT_ASSERT_MESSAGE("Test not implemented.", false); // :TODO: ADD MORE HERE

  PYLITH_METHOD_END;
} // testComputeLHSJacobianImplicit


// ----------------------------------------------------------------------
// Test computeLHSJacobianExplicit().
void
pylith::materials::TestIsotropicLinearElasticityPlaneStrain::testComputeLHSJacobianExplicit(void)
{ // testComputeLHSJacobianExplicit
  PYLITH_METHOD_BEGIN;

  CPPUNIT_ASSERT_MESSAGE("Test not implemented.", false); // :TODO: ADD MORE HERE

  PYLITH_METHOD_END;
} // testComputeLHSJacobianExplicit


// ----------------------------------------------------------------------
// Test updateStateVars().
void
pylith::materials::TestIsotropicLinearElasticityPlaneStrain::testUpdateStateVars(void)
{ // testUpdateStateVars
  PYLITH_METHOD_BEGIN;

  CPPUNIT_ASSERT_MESSAGE("Test not implemented.", false); // :TODO: ADD MORE HERE

  PYLITH_METHOD_END;  
} // testUpdateStateVars


// ----------------------------------------------------------------------
// Do minimal initilaization of test data.
void
pylith::materials::TestIsotropicLinearElasticityPlaneStrain::_initializeMin(void)
{ // _initializeMin
  PYLITH_METHOD_BEGIN;

  CPPUNIT_ASSERT(_material);
  CPPUNIT_ASSERT(_data);

  meshio::MeshIOAscii iohandler;
  CPPUNIT_ASSERT(_data->meshFilename);
  iohandler.filename(_data->meshFilename);
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

  _material->id(_data->materialId);
  _material->label(_data->materialLabel);
  _material->useInertia(_data->useInertia);
  _material->useBodyForce(_data->useBodyForce);
  _material->normalizer(normalizer);

  PYLITH_METHOD_END;
} // _initializeMin


// ----------------------------------------------------------------------
// Complete initilaization of test data.
void
pylith::materials::TestIsotropicLinearElasticityPlaneStrain::_initializeFull(void)
{ // _initializeFull
  PYLITH_METHOD_BEGIN;

  CPPUNIT_ASSERT(_material);
  CPPUNIT_ASSERT(_data);
  CPPUNIT_ASSERT(_mesh);

  _initializeMin();

  // Set auxiliary fields spatial database.
  delete _auxDB; _auxDB = new spatialdata::spatialdb::SimpleDB;CPPUNIT_ASSERT(_auxDB);
  spatialdata::spatialdb::SimpleIOAscii dbIO;
  CPPUNIT_ASSERT(_data->auxDBFilename);
  dbIO.filename(_data->auxDBFilename);
  _auxDB->ioHandler(&dbIO);
  _auxDB->label("IsotropicLinearElasciticityPlaneStrain auxiliary fields database");
  _auxDB->queryType(spatialdata::spatialdb::SimpleDB::LINEAR);
  _material->auxFieldsDB(_auxDB);

  // Create solution field.
  delete _solution1; _solution1 = new pylith::topology::Field(*_mesh);
  _solution1->label("solution1");
  CPPUNIT_ASSERT(_data->solnDiscretizations);
  const char* componentsDisp[2] = {"displacement_x", "displacement_y"};
  const char* componentsVel[2] = {"velocity_x", "velocity_y"};
  _solution1->subfieldAdd("displacement", componentsDisp, _data->dimension, topology::Field::VECTOR, _data->solnDiscretizations[0], _data->lengthScale);
  _solution1->subfieldAdd("velocity", componentsVel, _data->dimension, topology::Field::VECTOR, _data->solnDiscretizations[1], _data->lengthScale / _data->timeScale);
  _solution1->subfieldsSetup();
  _solution1->allocate();
  _solution1->zeroAll();

  delete _soln1DB; _soln1DB = new spatialdata::spatialdb::SimpleDB;CPPUNIT_ASSERT(_soln1DB);
  CPPUNIT_ASSERT(_data->soln1DBFilename);
  dbIO.filename(_data->soln1DBFilename);
  _soln1DB->ioHandler(&dbIO);
  _soln1DB->label("IsotropicLinearElasciticityPlaneStrain solution1 database");
  _soln1DB->queryType(spatialdata::spatialdb::SimpleDB::LINEAR);

  pylith::topology::FieldQuery querySoln1(*_solution1);
  querySoln1.queryFn("displacement", pylith::topology::FieldQuery::dbQueryGeneric);
  querySoln1.queryFn("velocity", pylith::topology::FieldQuery::dbQueryGeneric);
  querySoln1.openDB(_soln1DB, _data->lengthScale);
  querySoln1.queryDB();
  querySoln1.closeDB(_soln1DB);
  //_solution1->view("SOLUTION 1"); // DEBUGGING

  delete _solution2; _solution2 = new pylith::topology::Field(*_mesh);
  _solution2->cloneSection(*_solution1);
  _solution2->label("solution2");
  _solution2->allocate();
  _solution2->zeroAll();

  delete _soln2DB; _soln2DB = new spatialdata::spatialdb::SimpleDB;CPPUNIT_ASSERT(_soln2DB);
  CPPUNIT_ASSERT(_data->soln2DBFilename);
  dbIO.filename(_data->soln2DBFilename);
  _soln2DB->ioHandler(&dbIO);
  _soln2DB->label("IsotropicLinearElasciticityPlaneStrain solution2 database");
  _soln2DB->queryType(spatialdata::spatialdb::SimpleDB::LINEAR);

  pylith::topology::FieldQuery querySoln2(*_solution2);
  querySoln2.queryFn("displacement", pylith::topology::FieldQuery::dbQueryGeneric);
  querySoln2.queryFn("velocity", pylith::topology::FieldQuery::dbQueryGeneric);
  querySoln2.openDB(_soln2DB, _data->lengthScale);
  querySoln2.queryDB();
  querySoln2.closeDB(_soln2DB);
  //_solution2->view("SOLUTION 2"); // DEBUGGING
  
  delete _solution1Dot; _solution1Dot = new pylith::topology::Field(*_mesh);
  _solution1->label("solution1_dot");
  CPPUNIT_ASSERT(_data->solnDiscretizations);
  const char* componentsDispDot[2] = {"displacement_dot_x", "displacement_dot_y"};
  const char* componentsVelDot[2] = {"velocity_dot_x", "velocity_dot_y"};
  _solution1Dot->subfieldAdd("displacement_dot", componentsDispDot, _data->dimension, topology::Field::VECTOR, _data->solnDiscretizations[0], _data->lengthScale / _data->timeScale);
  _solution1Dot->subfieldAdd("velocity_dot", componentsVelDot, _data->dimension, topology::Field::VECTOR, _data->solnDiscretizations[1], _data->lengthScale / (_data->timeScale*_data->timeScale));
  _solution1Dot->subfieldsSetup();
  _solution1Dot->allocate();
  _solution1Dot->zeroAll();

  delete _solution2Dot; _solution2Dot = new pylith::topology::Field(*_mesh);
  _solution2Dot->cloneSection(*_solution1Dot);
  _solution2Dot->label("solution2");
  _solution2Dot->allocate();
  _solution2Dot->zeroAll();

  for (int i=0; i < _data->numAuxFields; ++i) {
    _material->auxFieldDiscretization(_data->auxFields[i], _data->auxDiscretizations[i]);
  } // for

  CPPUNIT_ASSERT(_solution1);
  _material->initialize(*_solution1);

  PYLITH_METHOD_END;
} // _initializeFull


// ----------------------------------------------------------------------
// Set field to zero on the boundary.
void
pylith::materials::TestIsotropicLinearElasticityPlaneStrain::_zeroBoundary(pylith::topology::Field* field)
{ // _zeroBoundary
  PYLITH_METHOD_BEGIN;

  CPPUNIT_ASSERT(field);
  CPPUNIT_ASSERT(_data);
  CPPUNIT_ASSERT(_data->boundaryLabel);

  PetscDM dmMesh = field->mesh().dmMesh();CPPUNIT_ASSERT(dmMesh);
  PetscDMLabel label = NULL;
  PetscIS pointIS = NULL;
  const PetscInt *points;
  PetscInt numPoints = 0;
  PetscBool hasLabel = PETSC_FALSE;
  PetscErrorCode err;
  err = DMHasLabel(dmMesh, _data->boundaryLabel, &hasLabel);CPPUNIT_ASSERT(!err);CPPUNIT_ASSERT(hasLabel);
  err = DMGetLabel(dmMesh, _data->boundaryLabel, &label);CPPUNIT_ASSERT(!err);
  err = DMLabelGetStratumIS(label, 1, &pointIS);CPPUNIT_ASSERT(!err);CPPUNIT_ASSERT(pointIS);
  err = ISGetLocalSize(pointIS, &numPoints);CPPUNIT_ASSERT(!err);
  err = ISGetIndices(pointIS, &points);CPPUNIT_ASSERT(!err);

  pylith::topology::VecVisitorMesh fieldVisitor(*field);
  PylithScalar* fieldArray = fieldVisitor.localArray();CPPUNIT_ASSERT(fieldArray);

  for (PetscInt p = 0; p < numPoints; ++p) {
    const PetscInt p_bc = points[p];

    const PylithInt off = fieldVisitor.sectionOffset(p_bc);
    const PylithInt dof = fieldVisitor.sectionDof(p_bc);
    for (PylithInt i=0; i < dof; ++i) {
      fieldArray[off+i] = 0.0;
    } // for
  } // for

  err = ISRestoreIndices(pointIS, &points);PYLITH_CHECK_ERROR(err);
  err = ISDestroy(&pointIS);PYLITH_CHECK_ERROR(err);

  PYLITH_METHOD_END;
} // _zeroBoundary




// End of file 
