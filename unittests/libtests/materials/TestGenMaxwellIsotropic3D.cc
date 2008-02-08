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
  CPPUNIT_ASSERT_EQUAL(&pylith::materials::GenMaxwellIsotropic3D::_updateStateElastic,
		       material._updateStateFn);

  material.useElasticBehavior(false);
  CPPUNIT_ASSERT_EQUAL(&pylith::materials::GenMaxwellIsotropic3D::_calcStressViscoelastic,
		       material._calcStressFn);
  CPPUNIT_ASSERT_EQUAL(&pylith::materials::GenMaxwellIsotropic3D::_calcElasticConstsViscoelastic,
		       material._calcElasticConstsFn);
  CPPUNIT_ASSERT_EQUAL(&pylith::materials::GenMaxwellIsotropic3D::_updateStateViscoelastic,
		       material._updateStateFn);
} // testUseElasticBehavior

// ----------------------------------------------------------------------
// Test usesUpdateState()
void
pylith::materials::TestGenMaxwellIsotropic3D::testUsesUpdateState(void)
{ // testUsesUpdateState
  GenMaxwellIsotropic3D material;
  CPPUNIT_ASSERT_EQUAL(true, material.usesUpdateState());
} // testUsesUpdateState

// ----------------------------------------------------------------------
// Test DBValues()
void
pylith::materials::TestGenMaxwellIsotropic3D::testDBValues(void)
{ // testDBValues
  GenMaxwellIsotropic3D material;
  GenMaxwellIsotropic3DElasticData data;
  material.useElasticBehavior(true);
  _testDBValues(&material, data);
} // testDBValues

// ----------------------------------------------------------------------
// Test parameters()
void
pylith::materials::TestGenMaxwellIsotropic3D::testParameters(void)
{ // testParameters
  GenMaxwellIsotropic3D material;
  GenMaxwellIsotropic3DElasticData data;
  material.useElasticBehavior(true);
  _testParameters(&material, data);
} // testParameters

// ----------------------------------------------------------------------
// Test _dbToParameters()
void
pylith::materials::TestGenMaxwellIsotropic3D::testDBToParameters(void)
{ // testDBToParameters
  GenMaxwellIsotropic3D material;
  GenMaxwellIsotropic3DElasticData data;
  material.useElasticBehavior(true);
  _testDBToParameters(&material, data);
} // testDBToParameters

// ----------------------------------------------------------------------
// Test calcDensity()
void
pylith::materials::TestGenMaxwellIsotropic3D::testCalcDensity(void)
{ // testCalcDensity
  GenMaxwellIsotropic3D material;
  GenMaxwellIsotropic3DElasticData data;
  material.useElasticBehavior(true);
  _testCalcDensity(&material, data);
} // testCalcDensity

// ----------------------------------------------------------------------
// Test calcStressElastic()
void
pylith::materials::TestGenMaxwellIsotropic3D::testCalcStressElastic(void)
{ // testCalcStressElastic
  GenMaxwellIsotropic3D material;
  GenMaxwellIsotropic3DElasticData data;
  material.useElasticBehavior(true);
  _testCalcStress(&material, data);
} // testCalcStressElastic

// ----------------------------------------------------------------------
// Test calcElasticConstsElastic()
void
pylith::materials::TestGenMaxwellIsotropic3D::testCalcElasticConstsElastic(void)
{ // testElasticConstsElastic
  GenMaxwellIsotropic3D material;
  GenMaxwellIsotropic3DElasticData data;
  material.useElasticBehavior(true);
  _testCalcElasticConsts(&material, data);
} // testElasticConstsElastic

// ----------------------------------------------------------------------
// Test updateStateElastic()
void
pylith::materials::TestGenMaxwellIsotropic3D::testUpdateStateElastic(void)
{ // testUpdateStateElastic
  GenMaxwellIsotropic3D material;
  GenMaxwellIsotropic3DElasticData data;

  const int_array& numParamValues = material._getNumParamValues();
  const int numParams = numParamValues.size();
  const int numParamsQuadPt = material._numParamsQuadPt;;

  const int tensorSize = 6;
  double_array totalStrain(tensorSize);
  for (int i=0; i < tensorSize; ++i)
    totalStrain[i] = i;

  const double meanStrain = 
    (totalStrain[0] + totalStrain[1] + totalStrain[2]) / 3.0;

  double_array parameters(numParamsQuadPt);
  double_array parametersE(numParamsQuadPt);
  for (int i=0, index=0; i < numParams; ++i)
    for (int j=0; j < numParamValues[i]; ++j, ++index) {
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
  
  material._updateState(&parameters[0], numParamsQuadPt, 
			&totalStrain[0], tensorSize);

  const double tolerance = 1.0e-06;
  for (int i=0; i < numParamsQuadPt; ++i)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(parametersE[i], parameters[i], tolerance);
} // testUpdateStateElastic

// ----------------------------------------------------------------------
// Test calcStressTimeDep()
void
pylith::materials::TestGenMaxwellIsotropic3D::testCalcStressTimeDep(void)
{ // testCalcStressTimeDep
  GenMaxwellIsotropic3D material;
  GenMaxwellIsotropic3DTimeDepData data;
  material.useElasticBehavior(false);
  double dt = 2.0e+5;
  material.timeStep(dt);
  _testCalcStress(&material, data);
} // testCalcStressTimeDep

// ----------------------------------------------------------------------
// Test calcElasticConstsTimeDep()
void
pylith::materials::TestGenMaxwellIsotropic3D::testCalcElasticConstsTimeDep(void)
{ // testElasticConstsTimeDep
  GenMaxwellIsotropic3D material;
  GenMaxwellIsotropic3DTimeDepData data;
  material.useElasticBehavior(false);
  double dt = 2.0e+5;
  material.timeStep(dt);
  _testCalcElasticConsts(&material, data);
} // testElasticConstsTimeDep

// ----------------------------------------------------------------------
// Test updateStateTimeDep()
void
pylith::materials::TestGenMaxwellIsotropic3D::testUpdateStateTimeDep(void)
{ // testUpdateStateTimeDep
  GenMaxwellIsotropic3D material;
  GenMaxwellIsotropic3DTimeDepData data;

  const int_array& numParamValues = material._getNumParamValues();
  const int numParams = numParamValues.size();
  const int numParamsQuadPt = material._numParamsQuadPt;;

  const int numMaxwellModels = 3;
  const int tensorSize = 6;
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
  for (int i=0; i < tensorSize; ++i) {
    totalStrainTpdt[i] = i;
    totalStrainT[i] = totalStrainTpdt[i] / 2.0;
    visStrainT[i] = totalStrainTpdt[i] / 4.0;
    visStrainT[i + tensorSize] = totalStrainTpdt[i] / 4.0;
    visStrainT[i + 2 * tensorSize] = totalStrainTpdt[i] / 4.0;
  } // for

  const double meanStrainTpdt = 
    (totalStrainTpdt[0] + totalStrainTpdt[1] + totalStrainTpdt[2]) / 3.0;
  const double meanStrainT = 
    (totalStrainT[0] + totalStrainT[1] + totalStrainT[2]) / 3.0;

  const double diag[] = { 1.0, 1.0, 1.0, 0.0, 0.0, 0.0 };

  double_array parameters(numParamsQuadPt);
  double_array parametersE(numParamsQuadPt);
  for (int i=0, index=0; i < numParams; ++i)
    for (int j=0; j < numParamValues[i]; ++j, ++index) {
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
  
  material._updateState(&parameters[0], numParamsQuadPt, 
			&totalStrainTpdt[0], tensorSize);

  const double tolerance = 1.0e-06;
  for (int i=0; i < numParamsQuadPt; ++i)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(parametersE[i], parameters[i], tolerance);
} // testUpdateStateTimeDep


// End of file 
