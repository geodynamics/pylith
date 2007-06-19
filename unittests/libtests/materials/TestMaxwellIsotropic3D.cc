// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
//                           Brad T. Aagaard
//                        U.S. Geological Survey
//
// {LicenseText}
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "TestMaxwellIsotropic3D.hh" // Implementation of class methods

#include "data/MaxwellIsotropic3DElasticData.hh" // USES MaxwellIsotropic3DElasticData
#include "data/MaxwellIsotropic3DTimeDepData.hh" // USES MaxwellIsotropic3DTimeDepData

#include "pylith/materials/MaxwellIsotropic3D.hh" // USES MaxwellIsotropic3D

#include <stdexcept> // TEMPORARY
// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::materials::TestMaxwellIsotropic3D );

// ----------------------------------------------------------------------
// Test timeStep()
void
pylith::materials::TestMaxwellIsotropic3D::testTimeStep(void)
{ // testTimeStep
  MaxwellIsotropic3D material;

  CPPUNIT_ASSERT_EQUAL(false, material._needNewJacobian);

  const double dt1 = 1.0;
  material.timeStep(dt1);
  CPPUNIT_ASSERT_EQUAL(dt1, material.Material::timeStep());
  CPPUNIT_ASSERT_EQUAL(false, material.needNewJacobian());

  const double dt2 = 2.0;
  material.timeStep(dt2);
  CPPUNIT_ASSERT_EQUAL(dt2, material.Material::timeStep());
  CPPUNIT_ASSERT_EQUAL(true, material.needNewJacobian());
} // testTimeStep

// ----------------------------------------------------------------------
// Test useElasticBehavior()
void
pylith::materials::TestMaxwellIsotropic3D::testUseElasticBehavior(void)
{ // testUseElasticBehavior
  MaxwellIsotropic3D material;

  material.useElasticBehavior(true);
  CPPUNIT_ASSERT_EQUAL(&pylith::materials::MaxwellIsotropic3D::_calcStressElastic,
		       material._calcStressFn);
  CPPUNIT_ASSERT_EQUAL(&pylith::materials::MaxwellIsotropic3D::_calcElasticConstsElastic,
		       material._calcElasticConstsFn);
  CPPUNIT_ASSERT_EQUAL(&pylith::materials::MaxwellIsotropic3D::_updateStateElastic,
		       material._updateStateFn);

  material.useElasticBehavior(false);
  CPPUNIT_ASSERT_EQUAL(&pylith::materials::MaxwellIsotropic3D::_calcStressViscoelastic,
		       material._calcStressFn);
  CPPUNIT_ASSERT_EQUAL(&pylith::materials::MaxwellIsotropic3D::_calcElasticConstsViscoelastic,
		       material._calcElasticConstsFn);
  CPPUNIT_ASSERT_EQUAL(&pylith::materials::MaxwellIsotropic3D::_updateStateViscoelastic,
		       material._updateStateFn);
} // testUseElasticBehavior

// ----------------------------------------------------------------------
// Test usesUpdateState()
void
pylith::materials::TestMaxwellIsotropic3D::testUsesUpdateState(void)
{ // testUsesUpdateState
  MaxwellIsotropic3D material;
  CPPUNIT_ASSERT_EQUAL(true, material.usesUpdateState());
} // testUsesUpdateState

// ----------------------------------------------------------------------
// Test DBValues()
void
pylith::materials::TestMaxwellIsotropic3D::testDBValues(void)
{ // testDBValues
  MaxwellIsotropic3D material;
  MaxwellIsotropic3DElasticData data;
  material.useElasticBehavior(true);
  _testDBValues(&material, data);
} // testDBValues

// ----------------------------------------------------------------------
// Test parameters()
void
pylith::materials::TestMaxwellIsotropic3D::testParameters(void)
{ // testParameters
  MaxwellIsotropic3D material;
  MaxwellIsotropic3DElasticData data;
  material.useElasticBehavior(true);
  _testParameters(&material, data);
} // testParameters

// ----------------------------------------------------------------------
// Test _dbToParameters()
void
pylith::materials::TestMaxwellIsotropic3D::testDBToParameters(void)
{ // testDBToParameters
  MaxwellIsotropic3D material;
  MaxwellIsotropic3DElasticData data;
  material.useElasticBehavior(true);
  _testDBToParameters(&material, data);
} // testDBToParameters

// ----------------------------------------------------------------------
// Test calcDensity()
void
pylith::materials::TestMaxwellIsotropic3D::testCalcDensity(void)
{ // testCalcDensity
  MaxwellIsotropic3D material;
  MaxwellIsotropic3DElasticData data;
  material.useElasticBehavior(true);
  _testCalcDensity(&material, data);
} // testCalcDensity

// ----------------------------------------------------------------------
// Test calcStressElastic()
void
pylith::materials::TestMaxwellIsotropic3D::testCalcStressElastic(void)
{ // testCalcStressElastic
  MaxwellIsotropic3D material;
  MaxwellIsotropic3DElasticData data;
  material.useElasticBehavior(true);
  _testCalcStress(&material, data);
} // testCalcStressElastic

// ----------------------------------------------------------------------
// Test calcElasticConstsElastic()
void
pylith::materials::TestMaxwellIsotropic3D::testCalcElasticConstsElastic(void)
{ // testElasticConstsElastic
  MaxwellIsotropic3D material;
  MaxwellIsotropic3DElasticData data;
  material.useElasticBehavior(true);
  _testCalcElasticConsts(&material, data);
} // testElasticConstsElastic

// ----------------------------------------------------------------------
// Test updateStateElastic()
void
pylith::materials::TestMaxwellIsotropic3D::testUpdateStateElastic(void)
{ // testUpdateStateElastic
  MaxwellIsotropic3D material;
  MaxwellIsotropic3DElasticData data;

  const int numParams = data.numParameters;
    
  const int tensorSize = 6;
  double_array totalStrain(tensorSize);
  for (int i=0; i < tensorSize; ++i)
    totalStrain[i] = i;

  const double meanStrain = (totalStrain[0] + totalStrain[1] + totalStrain[2])/3.0;

  std::vector<double_array> parameters(numParams);
  std::vector<double_array> paramdata(numParams);
  const int paramsSize [] = { 1, 1, 1, 1, 6, 6};
  for (int i=0; i < numParams; ++i) {
    parameters[i].resize(paramsSize[i]);
    paramdata[i].resize(paramsSize[i]);
    for (int j=0; j < paramsSize[i]; ++j)
      parameters[i][j] = i+j;
  } // for

  // Set up vector parameters, which are the only ones updated.
  for (int i=0; i< 3; ++i) {
    paramdata[4][i] = totalStrain[i];
    paramdata[5][i] = totalStrain[i] - meanStrain;
  } // for
  
  for (int i=3; i< 6; ++i) {
    paramdata[4][i] = totalStrain[i];
    paramdata[5][i] = totalStrain[i];
  } // for
  
  material._updateState(&parameters, totalStrain);

  const double tolerance = 1.0e-06;
  // Only test vector parameters
  for (int i=4; i < numParams; ++i)
    for (int j=0; j < paramsSize[i]; ++j)
      CPPUNIT_ASSERT_DOUBLES_EQUAL(paramdata[i][j], parameters[i][j], tolerance);
    
  for (int i=0; i < tensorSize; ++i)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(double(i), totalStrain[i], tolerance);
} // testUpdateStateElastic

// ----------------------------------------------------------------------
// Test calcStressTimeDep()
void
pylith::materials::TestMaxwellIsotropic3D::testCalcStressTimeDep(void)
{ // testCalcStressTimeDep
  MaxwellIsotropic3D material;
  MaxwellIsotropic3DTimeDepData data;
  material.useElasticBehavior(false);
  double dt = 2.0e5;
  material.timeStep(dt);
  _testCalcStress(&material, data);
} // testCalcStressTimeDep

// ----------------------------------------------------------------------
// Test calcElasticConstsTimeDep()
void
pylith::materials::TestMaxwellIsotropic3D::testCalcElasticConstsTimeDep(void)
{ // testElasticConstsTimeDep
  MaxwellIsotropic3D material;
  MaxwellIsotropic3DTimeDepData data;
  material.useElasticBehavior(false);
  double dt = 2.0e5;
  material.timeStep(dt);
  _testCalcElasticConsts(&material, data);
} // testElasticConstsTimeDep

// ----------------------------------------------------------------------
// Test updateStateTimeDep()
void
pylith::materials::TestMaxwellIsotropic3D::testUpdateStateTimeDep(void)
{ // testUpdateStateTimeDep
  MaxwellIsotropic3D material;
  MaxwellIsotropic3DTimeDepData data;

  // CHARLES: Setup parameters and totalStrain so this works.
#if 0
  const int numParams = data.numParameters;

  std::vector<double_array> parameters(numParams);
  const int paramsSize = 1;
  for (int i=0; i < numParams; ++i) {
    parameters[i].resize(numParams);
    for (int j=0; j < paramsSize; ++j)
      parameters[i][j] = i+j;
  } // for
    
  const int tensorSize = 9;
  double_array totalStrain(tensorSize);
  for (int i=0; i < tensorSize; ++i)
    totalStrain[i] = i;
  
  material._updateState(&parameters, totalStrain);

  const double tolerance = 1.0e-06;
  for (int i=0; i < numParams; ++i)
    for (int j=0; j < paramsSize; ++j)
      CPPUNIT_ASSERT_DOUBLES_EQUAL(double(i+j), parameters[i][j], tolerance);
    
  for (int i=0; i < tensorSize; ++i)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(double(i), totalStrain[i], tolerance);
#endif
  throw std::logic_error("Unit test not implemented.");
} // testUpdateStateTimeDep


// End of file 
