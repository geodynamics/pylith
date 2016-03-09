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
  _solutionDot = NULL;
  _db = NULL;
} // setUp


// ----------------------------------------------------------------------
// Deallocate testing data.
void
pylith::materials::TestIsotropicLinearElasticityPlaneStrain::tearDown(void)
{ // tearDown
  delete _solution; _solution = NULL;
  delete _solutionDot; _solutionDot = NULL;
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

  // Make sure DB query functions are set correctly.
  CPPUNIT_ASSERT_EQUAL(&pylith::topology::FieldQuery::dbQueryGeneric, _material->_auxFieldsQuery->queryFn("density"));
  CPPUNIT_ASSERT_EQUAL(&pylith::materials::Query::dbQueryMu2D, _material->_auxFieldsQuery->queryFn("mu"));
  CPPUNIT_ASSERT_EQUAL(&pylith::materials::Query::dbQueryLambda2D, _material->_auxFieldsQuery->queryFn("lambda"));
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

  CPPUNIT_ASSERT(_solution);
  CPPUNIT_ASSERT(_material);
  _material->_setFEKernelsRHSResidual(*_solution);

  PetscDS prob = NULL;
  PetscErrorCode err = DMGetDS(_solution->dmMesh(), &prob);CPPUNIT_ASSERT(!err);CPPUNIT_ASSERT(prob);

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

  CPPUNIT_ASSERT(_solution);
  CPPUNIT_ASSERT(_material);
  _material->_setFEKernelsRHSJacobian(*_solution);

  PetscDS prob = NULL;
  PetscErrorCode err = DMGetDS(_solution->dmMesh(), &prob);CPPUNIT_ASSERT(!err);CPPUNIT_ASSERT(prob);

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

  CPPUNIT_ASSERT(_solution);
  CPPUNIT_ASSERT(_material);
  _material->_setFEKernelsLHSResidual(*_solution);

  PetscDS prob = NULL;
  PetscErrorCode err = DMGetDS(_solution->dmMesh(), &prob);CPPUNIT_ASSERT(!err);CPPUNIT_ASSERT(prob);

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

  CPPUNIT_ASSERT(_solution);
  CPPUNIT_ASSERT(_material);
  _material->_setFEKernelsLHSJacobianImplicit(*_solution);

  PetscDS prob = NULL;
  PetscErrorCode err = DMGetDS(_solution->dmMesh(), &prob);CPPUNIT_ASSERT(!err);CPPUNIT_ASSERT(prob);

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

  CPPUNIT_ASSERT(_solution);
  CPPUNIT_ASSERT(_material);
  _material->_setFEKernelsLHSJacobianExplicit(*_solution);

  PetscDS prob = NULL;
  PetscErrorCode err = DMGetDS(_solution->dmMesh(), &prob);CPPUNIT_ASSERT(!err);CPPUNIT_ASSERT(prob);

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
  density.complete(); // Needed to populate global vector.

  //density.view("DENSITY"); // DEBUGGING

  // Check result
  CPPUNIT_ASSERT_EQUAL(std::string("density"), std::string(density.label()));
  CPPUNIT_ASSERT_EQUAL(_data->dimension, density.spaceDim());

  pylith::topology::FieldQuery queryDensity(density);
  queryDensity.queryFn("density", pylith::topology::FieldQuery::dbQueryGeneric);
  queryDensity.openDB(_db, _data->lengthScale);

  PylithReal norm = 0.0;
  PylithReal t = 0.0;
  const PetscDM dm = density.dmMesh();CPPUNIT_ASSERT(dm);
  PetscErrorCode err = DMPlexComputeL2Diff(dm, t, queryDensity.functions(), (void**)queryDensity.contextPtrs(), density.globalVector(), &norm);CPPUNIT_ASSERT(!err);
  queryDensity.closeDB(_db);

  const PylithReal tolerance = 1.0e-6;
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, norm, tolerance);


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
  query->openDB(_db, _data->lengthScale);

  PetscErrorCode err = DMPlexComputeL2Diff(dm, t, query->functions(), (void**)query->contextPtrs(), auxFields.globalVector(), &norm);CPPUNIT_ASSERT(!err);
  query->closeDB(_db);
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
  residualRHS.cloneSection(*_solution);
  residualRHS.label("residual RHS");
  residualRHS.allocate();
  residualRHS.zeroAll();

  pylith::topology::Field residualLHS(*_mesh);
  residualLHS.cloneSection(*_solution);
  residualLHS.label("residual LHS");
  residualLHS.allocate();
  residualLHS.zeroAll();

  pylith::topology::Field residual(*_mesh);
  residual.cloneSection(*_solution);
  residual.label("residual");
  residual.allocate();
  residual.zeroAll();

  pylith::topology::FieldQuery querySoln(*_solution);
  querySoln.queryFn("displacement", pylith::topology::FieldQuery::dbQueryGeneric);
  querySoln.queryFn("velocity", pylith::topology::FieldQuery::dbQueryGeneric);
  querySoln.openDB(_db, _data->lengthScale);
  querySoln.queryDB();
  querySoln.closeDB(_db);
  _solution->view("SOLUTION");

  pylith::topology::FieldQuery querySolnDot(*_solutionDot);
  querySolnDot.queryFn("displacement_dot", pylith::topology::FieldQuery::dbQueryGeneric);
  querySolnDot.queryFn("velocity_dot", pylith::topology::FieldQuery::dbQueryGeneric);
  querySolnDot.openDB(_db, _data->lengthScale);
  querySolnDot.queryDB();
  querySolnDot.closeDB(_db);
  _solutionDot->view("SOLUTION DOT");
  
  CPPUNIT_ASSERT(_material);
  PylithReal t = 1.0;
  PylithReal dt = 0.01;
  _material->computeRHSResidual(&residualRHS, t, dt, *_solution);
  _material->computeLHSResidual(&residualLHS, t, dt, *_solution, *_solutionDot);

  // Scatter local to global.
  residualRHS.complete();
  residualLHS.complete();

  residualRHS.view("RESIDUAL RHS");
  residualLHS.view("RESIDUAL LHS");

  PetscErrorCode err = VecWAXPY(residual.globalVector(), -1.0, residualLHS.globalVector(), residualRHS.globalVector());CPPUNIT_ASSERT(!err);

  PylithReal norm = 0.0;
  err = VecNorm(residual.globalVector(), NORM_2, &norm);CPPUNIT_ASSERT(!err);
  const PylithReal tolerance = 1.0e-6;
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, norm, tolerance);

  
  PYLITH_METHOD_END;
} // testComputeResidual


// ----------------------------------------------------------------------
// Test computeRHSJacobian().
void
pylith::materials::TestIsotropicLinearElasticityPlaneStrain::testComputeRHSJacobian(void)
{ // testComputeRHSJacobian
  PYLITH_METHOD_BEGIN;

  CPPUNIT_ASSERT_MESSAGE("Test not implemented.", false); // :TODO: ADD MORE HERE

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
  delete _db; _db = new spatialdata::spatialdb::SimpleDB;CPPUNIT_ASSERT(_db);
  spatialdata::spatialdb::SimpleIOAscii dbIO;
  dbIO.filename(_data->filenameAuxFieldsDB);
  _db->ioHandler(&dbIO);
  _db->label("IsotropicLinearElasciticityPlaneStrain auxiliary fields");
  _db->queryType(spatialdata::spatialdb::SimpleDB::LINEAR);
  _material->auxFieldsDB(_db);

  // Create solution field.
  delete _solution; _solution = new pylith::topology::Field(*_mesh);
  _solution->label("solution");
  CPPUNIT_ASSERT(_data->solnDiscretizations);
  const char* componentsDisp[2] = {"displacement_x", "displacement_y"};
  const char* componentsVel[2] = {"velocity_x", "velocity_y"};
  _solution->subfieldAdd("displacement", componentsDisp, _data->dimension, topology::Field::VECTOR, _data->solnDiscretizations[0], _data->lengthScale);
  _solution->subfieldAdd("velocity", componentsVel, _data->dimension, topology::Field::VECTOR, _data->solnDiscretizations[1], _data->lengthScale / _data->timeScale);
  _solution->subfieldsSetup();
  _solution->allocate();
  _solution->zeroAll();

  delete _solutionDot; _solutionDot = new pylith::topology::Field(*_mesh);
  _solution->label("solution_dot");
  CPPUNIT_ASSERT(_data->solnDiscretizations);
  const char* componentsDispDot[2] = {"displacement_dot_x", "displacement_dot_y"};
  const char* componentsVelDot[2] = {"velocity_dot_x", "velocity_dot_y"};
  _solutionDot->subfieldAdd("displacement_dot", componentsDispDot, _data->dimension, topology::Field::VECTOR, _data->solnDiscretizations[0], _data->lengthScale / _data->timeScale);
  _solutionDot->subfieldAdd("velocity_dot", componentsVelDot, _data->dimension, topology::Field::VECTOR, _data->solnDiscretizations[1], _data->lengthScale / (_data->timeScale*_data->timeScale));
  _solutionDot->subfieldsSetup();
  _solutionDot->allocate();
  _solutionDot->zeroAll();

  for (int i=0; i < _data->numAuxFields; ++i) {
    _material->auxFieldDiscretization(_data->auxFields[i], _data->auxDiscretizations[i]);
  } // for

  CPPUNIT_ASSERT(_solution);
  _material->initialize(*_solution);

  PYLITH_METHOD_END;
} // _initializeFull


// End of file 
