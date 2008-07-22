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

#include "TestGenMaxwellIsotropic3D.hh" // Implementation of class methods

#include "data/GenMaxwellIsotropic3DElasticData.hh" // USES GenMaxwellIsotropic3DElasticData
#include "data/GenMaxwellIsotropic3DTimeDepData.hh" // USES GenMaxwellIsotropic3DTimeDepData

#include "pylith/materials/GenMaxwellIsotropic3D.hh" // USES GenMaxwellIsotropic3D

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::materials::TestGenMaxwellIsotropic3D );

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::materials::TestGenMaxwellIsotropic3D::setUp(void)
{ // setUp
  _material = new GenMaxwellIsotropic3D();
  _matElastic = new GenMaxwellIsotropic3D();
  _data = new GenMaxwellIsotropic3DElasticData();
  _dataElastic = new GenMaxwellIsotropic3DElasticData();
} // setUp

// ----------------------------------------------------------------------
// Test timeStep()
void
pylith::materials::TestGenMaxwellIsotropic3D::testTimeStep(void)
{ // testTimeStep
  GenMaxwellIsotropic3D material;

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
pylith::materials::TestGenMaxwellIsotropic3D::testUseElasticBehavior(void)
{ // testUseElasticBehavior
  GenMaxwellIsotropic3D material;

  material.useElasticBehavior(true);
  CPPUNIT_ASSERT_EQUAL(&pylith::materials::GenMaxwellIsotropic3D::_calcStressElastic,
		       material._calcStressFn);
  CPPUNIT_ASSERT_EQUAL(&pylith::materials::GenMaxwellIsotropic3D::_calcElasticConstsElastic,
		       material._calcElasticConstsFn);
  CPPUNIT_ASSERT_EQUAL(&pylith::materials::GenMaxwellIsotropic3D::_updatePropertiesElastic,
		       material._updatePropertiesFn);

  material.useElasticBehavior(false);
  CPPUNIT_ASSERT_EQUAL(&pylith::materials::GenMaxwellIsotropic3D::_calcStressViscoelastic,
		       material._calcStressFn);
  CPPUNIT_ASSERT_EQUAL(&pylith::materials::GenMaxwellIsotropic3D::_calcElasticConstsViscoelastic,
		       material._calcElasticConstsFn);
  CPPUNIT_ASSERT_EQUAL(&pylith::materials::GenMaxwellIsotropic3D::_updatePropertiesViscoelastic,
		       material._updatePropertiesFn);
} // testUseElasticBehavior

// ----------------------------------------------------------------------
// Test usesUpdateProperties()
void
pylith::materials::TestGenMaxwellIsotropic3D::testUsesUpdateProperties(void)
{ // testUsesUpdateProperties
  GenMaxwellIsotropic3D material;
  CPPUNIT_ASSERT_EQUAL(true, material.usesUpdateProperties());
} // testUsesUpdateProperties

// ----------------------------------------------------------------------
// Test calcStressElastic()
void
pylith::materials::TestGenMaxwellIsotropic3D::testCalcStressElastic(void)
{ // testCalcStressElastic
  CPPUNIT_ASSERT(0 != _matElastic);
  _matElastic->useElasticBehavior(true);

  test_calcStress();
} // testCalcStressElastic

// ----------------------------------------------------------------------
// Test calcElasticConstsElastic()
void
pylith::materials::TestGenMaxwellIsotropic3D::testCalcElasticConstsElastic(void)
{ // testElasticConstsElastic
  CPPUNIT_ASSERT(0 != _matElastic);
  _matElastic->useElasticBehavior(true);

  test_calcElasticConsts();
} // testElasticConstsElastic

// ----------------------------------------------------------------------
// Test updatePropertiesElastic()
void
pylith::materials::TestGenMaxwellIsotropic3D::testUpdatePropertiesElastic(void)
{ // testUpdatePropertiesElastic
  GenMaxwellIsotropic3D material;
  GenMaxwellIsotropic3DElasticData data;

  const int numParams = data.numParameters;
  const int numParamsQuadPt = data.numParamsQuadPt;

  const int tensorSize = 6;
  const int initialStateSize = 6;
  double_array totalStrain(tensorSize);
  double_array initialState(initialStateSize);
  for (int i=0; i < tensorSize; ++i) {
    totalStrain[i] = i;
    initialState[i] = 0;
  } // for

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
  const int numMaxwellModels = 3;
  const int pidStrainT = 3 + 2 * numMaxwellModels;
  const int pidVisStrain = pidStrainT + tensorSize;
  for (int i=0; i < 3; ++i) {
    parametersE[pidStrainT + i] = totalStrain[i];
    for (int j=0; j < numMaxwellModels; ++j)
      parametersE[pidVisStrain + i + j * tensorSize] =
	totalStrain[i] - meanStrain;
  } // for
  
  for (int i=3; i < 6; ++i) {
    parametersE[pidStrainT + i] = totalStrain[i];
    for (int j=0; j < numMaxwellModels; ++j)
      parametersE[pidVisStrain + i + j * tensorSize] =
	totalStrain[i];
  } // for
  
  material._updateProperties(&parameters[0], numParamsQuadPt, 
			&totalStrain[0], tensorSize,
			&initialState[0], initialStateSize);

  const double tolerance = 1.0e-06;
  for (int i=0; i < numParamsQuadPt; ++i)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(parametersE[i], parameters[i], tolerance);
} // testUpdatePropertiesElastic

// ----------------------------------------------------------------------
// Test calcStressTimeDep()
void
pylith::materials::TestGenMaxwellIsotropic3D::testCalcStressTimeDep(void)
{ // testCalcStressTimeDep
  CPPUNIT_ASSERT(0 != _matElastic);
  _matElastic->useElasticBehavior(false);

  delete _dataElastic; _dataElastic = new GenMaxwellIsotropic3DTimeDepData();

  double dt = 2.0e+5;
  _matElastic->timeStep(dt);
  test_calcStress();
} // testCalcStressTimeDep

// ----------------------------------------------------------------------
// Test calcElasticConstsTimeDep()
void
pylith::materials::TestGenMaxwellIsotropic3D::testCalcElasticConstsTimeDep(void)
{ // testElasticConstsTimeDep
  CPPUNIT_ASSERT(0 != _matElastic);
  _matElastic->useElasticBehavior(false);

  delete _dataElastic; _dataElastic = new GenMaxwellIsotropic3DTimeDepData();

  double dt = 2.0e+5;
  _matElastic->timeStep(dt);
  test_calcElasticConsts();
} // testElasticConstsTimeDep

// ----------------------------------------------------------------------
// Test updatePropertiesTimeDep()
void
pylith::materials::TestGenMaxwellIsotropic3D::testUpdatePropertiesTimeDep(void)
{ // testUpdatePropertiesTimeDep
  GenMaxwellIsotropic3D material;
  GenMaxwellIsotropic3DTimeDepData data;

  const int numParams = data.numParameters;
  const int numParamsQuadPt = data.numParamsQuadPt;

  const int numMaxwellModels = 3;
  const int tensorSize = 6;
  const int initialStateSize = 6;
  const double mu = 3.0e10;

  material.useElasticBehavior(false);
  const double dt = 2.0e+5;
  material.timeStep(dt);
  const double shearRatio[] = {0.2, 0.3, 0.4};
  const double viscosity[] = {1.0e+18, 2.0e17, 3.0e19};
  double maxwellTime[] = {0.0, 0.0, 0.0};
  for (int model = 0; model < numMaxwellModels; ++model)
    maxwellTime[model] = viscosity[model]/(mu * shearRatio[model]);
    
  double_array totalStrainTpdt(tensorSize);
  double_array totalStrainT(tensorSize);
  double_array visStrainT(numMaxwellModels * tensorSize);
  double_array initialState(initialStateSize);
  for (int i=0; i < tensorSize; ++i) {
    totalStrainTpdt[i] = i;
    totalStrainT[i] = totalStrainTpdt[i] / 2.0;
    visStrainT[i] = totalStrainTpdt[i] / 4.0;
    visStrainT[i + tensorSize] = totalStrainTpdt[i] / 4.0;
    visStrainT[i + 2 * tensorSize] = totalStrainTpdt[i] / 4.0;
    initialState[i] = 0;
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


  const int pidMuTot = 1;
  const int pidShearRatio = 3;
  const int pidMaxwellTime = pidShearRatio + numMaxwellModels;
  const int pidStrainT = pidMaxwellTime + numMaxwellModels;
  const int pidVisStrain = pidStrainT + tensorSize;

  parameters[pidMuTot] = mu;
  parametersE[pidMuTot] = mu;

  double dq[] = {0.0, 0.0, 0.0};
  for (int model = 0; model < numMaxwellModels; ++model) {
    parameters[pidShearRatio + model] = shearRatio[model];
    parameters[pidMaxwellTime + model] = maxwellTime[model];
    parametersE[pidShearRatio + model] = shearRatio[model];
    parametersE[pidMaxwellTime + model] = maxwellTime[model];
    dq[model] = maxwellTime[model] * (1.0 - exp(-dt/maxwellTime[model]))/dt;
  }

  double devStrainTpdt = 0.0;
  double devStrainT = 0.0;
  double deltaStrain = 0.0;
  double visStrain = 0.0;

  for (int iComp = 0; iComp < tensorSize; ++iComp) {
    devStrainTpdt = totalStrainTpdt[iComp] - diag[iComp]*meanStrainTpdt;
    devStrainT = totalStrainT[iComp] - diag[iComp]*meanStrainT;
    deltaStrain = devStrainTpdt - devStrainT;
    parameters[pidStrainT + iComp] = totalStrainT[iComp];
    parametersE[pidStrainT + iComp] = totalStrainTpdt[iComp];
    for (int model = 0; model < numMaxwellModels; ++model) {
      int index = iComp + model * tensorSize;
      parametersE[pidVisStrain + index] =
	exp(-dt/maxwellTime[model]) *
	visStrainT[index] + dq[model] * deltaStrain;
      parameters[pidVisStrain + index] = visStrainT[index];
    } // for
  } // for
  
  material._updateProperties(&parameters[0], numParamsQuadPt, 
			&totalStrainTpdt[0], tensorSize,
			&initialState[0], initialStateSize);

  const double tolerance = 1.0e-06;
  for (int i=0; i < numParamsQuadPt; ++i)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(parametersE[i], parameters[i], tolerance);
} // testUpdatePropertiesTimeDep

// ----------------------------------------------------------------------
// Test _stableTimeStepImplicit()
void
pylith::materials::TestGenMaxwellIsotropic3D::test_stableTimeStepImplicit(void)
{ // test_stableTimeStepImplicit
  CPPUNIT_ASSERT(0 != _matElastic);

  delete _dataElastic; _dataElastic = new GenMaxwellIsotropic3DTimeDepData();

  TestElasticMaterial::test_stableTimeStepImplicit();
} // test_stableTimeStepImplicit


// End of file 
