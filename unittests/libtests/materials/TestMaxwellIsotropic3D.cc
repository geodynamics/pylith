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

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::materials::TestMaxwellIsotropic3D );

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::materials::TestMaxwellIsotropic3D::setUp(void)
{ // setUp
  _material = new MaxwellIsotropic3D();
  _matElastic = new MaxwellIsotropic3D();
  _data = new MaxwellIsotropic3DElasticData();
  _dataElastic = new MaxwellIsotropic3DElasticData();
} // setUp

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
  CPPUNIT_ASSERT_EQUAL(&pylith::materials::MaxwellIsotropic3D::_updatePropertiesElastic,
		       material._updatePropertiesFn);

  material.useElasticBehavior(false);
  CPPUNIT_ASSERT_EQUAL(&pylith::materials::MaxwellIsotropic3D::_calcStressViscoelastic,
		       material._calcStressFn);
  CPPUNIT_ASSERT_EQUAL(&pylith::materials::MaxwellIsotropic3D::_calcElasticConstsViscoelastic,
		       material._calcElasticConstsFn);
  CPPUNIT_ASSERT_EQUAL(&pylith::materials::MaxwellIsotropic3D::_updatePropertiesViscoelastic,
		       material._updatePropertiesFn);
} // testUseElasticBehavior

// ----------------------------------------------------------------------
// Test usesUpdateProperties()
void
pylith::materials::TestMaxwellIsotropic3D::testUsesUpdateProperties(void)
{ // testUsesUpdateProperties
  MaxwellIsotropic3D material;
  CPPUNIT_ASSERT_EQUAL(true, material.usesUpdateProperties());
} // testUsesUpdateProperties

// ----------------------------------------------------------------------
// Test calcStressElastic()
void
pylith::materials::TestMaxwellIsotropic3D::testCalcStressElastic(void)
{ // testCalcStressElastic
  CPPUNIT_ASSERT(0 != _matElastic);
  _matElastic->useElasticBehavior(true);

  test_calcStress();
} // testCalcStressElastic

// ----------------------------------------------------------------------
// Test calcElasticConstsElastic()
void
pylith::materials::TestMaxwellIsotropic3D::testCalcElasticConstsElastic(void)
{ // testElasticConstsElastic
  CPPUNIT_ASSERT(0 != _matElastic);
  _matElastic->useElasticBehavior(true);

  test_calcElasticConsts();
} // testElasticConstsElastic

// ----------------------------------------------------------------------
// Test updatePropertiesElastic()
void
pylith::materials::TestMaxwellIsotropic3D::testUpdatePropertiesElastic(void)
{ // testUpdatePropertiesElastic
  MaxwellIsotropic3D material;
  MaxwellIsotropic3DElasticData data;

  const int numParams = data.numParameters;
  const int numParamsQuadPt = data.numParamsQuadPt;

  const int tensorSize = 6;
  double_array totalStrain(tensorSize);
  for (int i=0; i < tensorSize; ++i)
    totalStrain[i] = i;

  const double meanStrain = 
    (totalStrain[0] + totalStrain[1] + totalStrain[2]) / 3.0;

  double_array parameters(numParamsQuadPt);
  double_array parametersE(numParamsQuadPt);
  for (int i=0, index=0; i < numParams; ++i)
    for (int j=0; j < data.numParamValues[i]; ++j, ++index) {
      parametersE[index] = i+j;
      parameters[index] = i+j;
    } // for

  // Set up vector parameters, which are the only ones updated.
  const int pidStrainT = 4;
  const int pidVisStrain = pidStrainT + 6;
  for (int i=0; i < 3; ++i) {
    parametersE[pidStrainT+i] = totalStrain[i];
    parametersE[pidVisStrain+i] = totalStrain[i] - meanStrain;
  } // for
  
  for (int i=3; i < 6; ++i) {
    parametersE[pidStrainT+i] = totalStrain[i];
    parametersE[pidVisStrain+i] = totalStrain[i];
  } // for
  
  material._updateProperties(&parameters[0], numParamsQuadPt, 
			&totalStrain[0], tensorSize);

  const double tolerance = 1.0e-06;
  for (int i=0; i < numParamsQuadPt; ++i)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(parametersE[i], parameters[i], tolerance);
} // testUpdatePropertiesElastic

// ----------------------------------------------------------------------
// Test calcStressTimeDep()
void
pylith::materials::TestMaxwellIsotropic3D::testCalcStressTimeDep(void)
{ // testCalcStressTimeDep
  CPPUNIT_ASSERT(0 != _matElastic);
  _matElastic->useElasticBehavior(false);

  delete _dataElastic; _dataElastic = new MaxwellIsotropic3DTimeDepData();

  double dt = 2.0e+5;
  _matElastic->timeStep(dt);
  test_calcStress();
} // testCalcStressTimeDep

// ----------------------------------------------------------------------
// Test calcElasticConstsTimeDep()
void
pylith::materials::TestMaxwellIsotropic3D::testCalcElasticConstsTimeDep(void)
{ // testElasticConstsTimeDep
  CPPUNIT_ASSERT(0 != _matElastic);
  _matElastic->useElasticBehavior(false);

  delete _dataElastic; _dataElastic = new MaxwellIsotropic3DTimeDepData();

  double dt = 2.0e+5;
  _matElastic->timeStep(dt);
  test_calcElasticConsts();
} // testElasticConstsTimeDep

// ----------------------------------------------------------------------
// Test updatePropertiesTimeDep()
void
pylith::materials::TestMaxwellIsotropic3D::testUpdatePropertiesTimeDep(void)
{ // testUpdatePropertiesTimeDep
  MaxwellIsotropic3D material;
  MaxwellIsotropic3DTimeDepData data;

  const int numParams = data.numParameters;
  const int numParamsQuadPt = data.numParamsQuadPt;

  material.useElasticBehavior(false);
  const double dt = 2.0e+5;
  material.timeStep(dt);
  const double viscosity = 1.0e+18;
  const double mu = 3.0e+10;
  const double maxwelltime = viscosity / mu;
    
  const int tensorSize = 6;
  double_array totalStrainTpdt(tensorSize);
  double_array totalStrainT(tensorSize);
  double_array visStrainT(tensorSize);
  for (int i=0; i < tensorSize; ++i) {
    totalStrainTpdt[i] = i;
    totalStrainT[i] = totalStrainTpdt[i] / 2.0;
    visStrainT[i] = totalStrainTpdt[i] / 4.0;
  } // for

  const double meanStrainTpdt = 
    (totalStrainTpdt[0] + totalStrainTpdt[1] + totalStrainTpdt[2]) / 3.0;
  const double meanStrainT = 
    (totalStrainT[0] + totalStrainT[1] + totalStrainT[2]) / 3.0;

  const double diag[] = { 1.0, 1.0, 1.0, 0.0, 0.0, 0.0 };

  double_array parameters(numParamsQuadPt);
  double_array parametersE(numParamsQuadPt);
  for (int i=0, index=0; i < numParams; ++i)
    for (int j=0; j < data.numParamValues[i]; ++j, ++index) {
      parametersE[index] = i+j;
      parameters[index] = i+j;
    } // for

  const int pidMaxwellTime = 3;
  const int pidStrainT = pidMaxwellTime + 1;
  const int pidVisStrain = pidStrainT + 6;

  parameters[pidMaxwellTime] = maxwelltime;
  parametersE[pidMaxwellTime] = maxwelltime;

  const double dq = maxwelltime*(1.0-exp(-dt/maxwelltime))/dt;
  const double expFac = exp(-dt/maxwelltime);
  double devStrainTpdt = 0.0;
  double devStrainT = 0.0;

  for (int i=0; i < tensorSize; ++i) {
    devStrainTpdt = totalStrainTpdt[i] - diag[i]*meanStrainTpdt;
    devStrainT = totalStrainT[i] - diag[i]*meanStrainT;
    parameters[pidStrainT+i] = totalStrainT[i];
    parameters[pidVisStrain+i] = visStrainT[i];
    parametersE[pidStrainT+i] = totalStrainTpdt[i];
    parametersE[pidVisStrain+i] = 
      expFac * visStrainT[i] + dq * (devStrainTpdt - devStrainT);
  } //for
  
  material._updateProperties(&parameters[0], numParamsQuadPt, 
			&totalStrainTpdt[0], tensorSize);

  const double tolerance = 1.0e-06;
  for (int i=0; i < numParamsQuadPt; ++i)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(parametersE[i], parameters[i], tolerance);
} // testUpdatePropertiesTimeDep

// ----------------------------------------------------------------------
// Test _stableTimeStepImplicit()
void
pylith::materials::TestMaxwellIsotropic3D::test_stableTimeStepImplicit(void)
{ // test_stableTimeStepImplicit
  CPPUNIT_ASSERT(0 != _matElastic);

  delete _dataElastic; _dataElastic = new MaxwellIsotropic3DTimeDepData();

  TestElasticMaterial::test_stableTimeStepImplicit();
} // test_stableTimeStepImplicit


// End of file 
