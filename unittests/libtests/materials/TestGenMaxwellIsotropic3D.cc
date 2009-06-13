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
  setupNormalizer();
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
  CPPUNIT_ASSERT_EQUAL(&pylith::materials::GenMaxwellIsotropic3D::_updateStateVarsElastic,
		       material._updateStateVarsFn);

  material.useElasticBehavior(false);
  CPPUNIT_ASSERT_EQUAL(&pylith::materials::GenMaxwellIsotropic3D::_calcStressViscoelastic,
		       material._calcStressFn);
  CPPUNIT_ASSERT_EQUAL(&pylith::materials::GenMaxwellIsotropic3D::_calcElasticConstsViscoelastic,
		       material._calcElasticConstsFn);
  CPPUNIT_ASSERT_EQUAL(&pylith::materials::GenMaxwellIsotropic3D::_updateStateVarsViscoelastic,
		       material._updateStateVarsFn);
} // testUseElasticBehavior

// ----------------------------------------------------------------------
// Test usesHasStateVars()
void
pylith::materials::TestGenMaxwellIsotropic3D::testHasStateVars(void)
{ // testHasStateVars
  GenMaxwellIsotropic3D material;
  CPPUNIT_ASSERT_EQUAL(true, material.hasStateVars());
} // testHasStateVars

// ----------------------------------------------------------------------
// Test calcStressElastic()
void
pylith::materials::TestGenMaxwellIsotropic3D::test_calcStressElastic(void)
{ // testCalcStressElastic
  CPPUNIT_ASSERT(0 != _matElastic);
  _matElastic->useElasticBehavior(true);

  test_calcStress();
} // testCalcStressElastic

// ----------------------------------------------------------------------
// Test calcElasticConstsElastic()
void
pylith::materials::TestGenMaxwellIsotropic3D::test_calcElasticConstsElastic(void)
{ // testElasticConstsElastic
  CPPUNIT_ASSERT(0 != _matElastic);
  _matElastic->useElasticBehavior(true);

  test_calcElasticConsts();
} // testElasticConstsElastic

// ----------------------------------------------------------------------
// Test updateStateVarsElastic()
void
pylith::materials::TestGenMaxwellIsotropic3D::test_updateStateVarsElastic(void)
{ // testUpdateStateVarsElastic
  // :TODO: Use TestElasticMaterial::test_updateStateVars
  // instead. This requires moving the calculation of the expected
  // state vars below to the Python code (where it belongs) and
  // setting the stateVarsUpdate attribute in the Python object.

  GenMaxwellIsotropic3D material;
  material.useElasticBehavior(true);

  const bool computeStateVars = true;
  
  const int numLocs = _dataElastic->numLocs;
  const int numPropsQuadPt = _dataElastic->numPropsQuadPt;
  const int numVarsQuadPt = _dataElastic->numVarsQuadPt;
  const int tensorSize = material.tensorSize();
  
  double_array stress(tensorSize);
  double_array properties(numPropsQuadPt);
  double_array stateVars(numVarsQuadPt);
  double_array strain(tensorSize);
  double_array initialStress(tensorSize);
  double_array initialStrain(tensorSize);
  
  for (int iLoc=0; iLoc < numLocs; ++iLoc) {
    memcpy(&properties[0], &_dataElastic->properties[iLoc*numPropsQuadPt],
	   numPropsQuadPt*sizeof(double));
    memcpy(&stateVars[0], &_dataElastic->stateVars[iLoc*numVarsQuadPt],
	   numVarsQuadPt*sizeof(double));
    memcpy(&strain[0], &_dataElastic->strain[iLoc*tensorSize],
	   tensorSize*sizeof(double));
    memcpy(&initialStress[0], &_dataElastic->initialStress[iLoc*tensorSize],
	   tensorSize*sizeof(double));
    memcpy(&initialStrain[0], &_dataElastic->initialStrain[iLoc*tensorSize],
	   tensorSize*sizeof(double));

    const double meanStrain = (strain[0] + strain[1] + strain[2]) / 3.0;
    
    // Compute expected state variables
    double_array stateVarsE(numVarsQuadPt);
    const int s_totalStrain = 0;
    const int s_viscousStrain = s_totalStrain + tensorSize;

    // State variable 'total_strain' should match 'strain'
    for (int i=0; i < tensorSize; ++i) 
      stateVarsE[s_totalStrain+i] = strain[i];
    
    // State variable 'viscous_strain'
    const double diag[] = { 1.0, 1.0, 1.0, 0.0, 0.0, 0.0 };
    const int numMaxwellModels = 3;
    for (int imodel=0; imodel < numMaxwellModels; ++imodel)
      for (int i=0; i < tensorSize; ++i)
	stateVarsE[s_viscousStrain+imodel*tensorSize+i] =
	  strain[i] - diag[i]*meanStrain;

    material._updateStateVars(&stateVars[0], stateVars.size(), 
			      &properties[0], properties.size(),
			      &strain[0], strain.size(),
			      &initialStress[0], initialStress.size(),
			      &initialStrain[0], initialStrain.size());

    const double tolerance = 1.0e-06;
    for (int i=0; i < numVarsQuadPt; ++i)
      CPPUNIT_ASSERT_DOUBLES_EQUAL(stateVarsE[i], stateVars[i], tolerance);
  } // for
} // testUpdateStateVarsElastic

// ----------------------------------------------------------------------
// Test calcStressTimeDep()
void
pylith::materials::TestGenMaxwellIsotropic3D::test_calcStressTimeDep(void)
{ // testCalcStressTimeDep
#if 0
  CPPUNIT_ASSERT(0 != _matElastic);
  _matElastic->useElasticBehavior(false);

  delete _dataElastic; _dataElastic = new GenMaxwellIsotropic3DTimeDepData();

  double dt = 2.0e+5;
  _matElastic->timeStep(dt);
  test_calcStress();
#endif
} // testCalcStressTimeDep

// ----------------------------------------------------------------------
// Test calcElasticConstsTimeDep()
void
pylith::materials::TestGenMaxwellIsotropic3D::test_calcElasticConstsTimeDep(void)
{ // testElasticConstsTimeDep
#if 0
  CPPUNIT_ASSERT(0 != _matElastic);
  _matElastic->useElasticBehavior(false);

  delete _dataElastic; _dataElastic = new GenMaxwellIsotropic3DTimeDepData();

  double dt = 2.0e+5;
  _matElastic->timeStep(dt);
  test_calcElasticConsts();
#endif
} // testElasticConstsTimeDep

// ----------------------------------------------------------------------
// Test updateStateVarsTimeDep()
void
pylith::materials::TestGenMaxwellIsotropic3D::test_updateStateVarsTimeDep(void)
{ // testUpdateStateVarsTimeDep
#if 0
  // :TODO: Use TestElasticMaterial::test_updateStateVars
  // instead. This requires moving the calculation of the expected
  // state vars below to the Python code (where it belongs) and
  // setting the stateVarsUpdate attribute in the Python object.

  MaxwellIsotropic3D material;
  material.useElasticBehavior(false);

  delete _dataElastic; _dataElastic = new MaxwellIsotropic3DTimeDepData();

  const double dt = 2.0e+5;
  material.timeStep(dt);

  const bool computeStateVars = true;
  
  const int numLocs = _dataElastic->numLocs;
  const int numPropsQuadPt = _dataElastic->numPropsQuadPt;
  const int numVarsQuadPt = _dataElastic->numVarsQuadPt;
  const int tensorSize = material.tensorSize();
  
  double_array stress(tensorSize);
  double_array properties(numPropsQuadPt);
  double_array stateVars(numVarsQuadPt);
  double_array strain(tensorSize);
  double_array initialStress(tensorSize);
  double_array initialStrain(tensorSize);
  
  for (int iLoc=0; iLoc < numLocs; ++iLoc) {
    memcpy(&properties[0], &_dataElastic->properties[iLoc*numPropsQuadPt],
	   numPropsQuadPt*sizeof(double));
    memcpy(&stateVars[0], &_dataElastic->stateVars[iLoc*numVarsQuadPt],
	   numVarsQuadPt*sizeof(double));
    memcpy(&strain[0], &_dataElastic->strain[iLoc*tensorSize],
	   tensorSize*sizeof(double));
    memcpy(&initialStress[0], &_dataElastic->initialStress[iLoc*tensorSize],
	   tensorSize*sizeof(double));
    memcpy(&initialStrain[0], &_dataElastic->initialStrain[iLoc*tensorSize],
	   tensorSize*sizeof(double));

    const double meanStrain = (strain[0] + strain[1] + strain[2]) / 3.0;
    
    // Compute expected state variables
    double_array stateVarsE(numVarsQuadPt);
    const int s_totalStrain = 0;
    const int s_viscousStrain = s_totalStrain + tensorSize;

    // State variable 'total_strain' should match 'strain'
    for (int i=0; i < tensorSize; ++i) 
      stateVarsE[s_totalStrain+i] = strain[i];
    
    // State variable 'viscous_strain'
    const double meanStrainTpdt = 
      (strain[0] + strain[1] + strain[2]) / 3.0;
    const double meanStrainT = 
      (stateVars[s_totalStrain+0] + 
       stateVars[s_totalStrain+1] + 
       stateVars[s_totalStrain+2]) / 3.0;

    const int p_maxwellTime = 3;
    const double maxwellTime = properties[p_maxwellTime];
    const double dq = maxwellTime*(1.0-exp(-dt/maxwellTime))/dt;
    const double expFac = exp(-dt/maxwellTime);
    double devStrainTpdt = 0.0;
    double devStrainT = 0.0;

    const double diag[] = { 1.0, 1.0, 1.0, 0.0, 0.0, 0.0 };
    for (int i=0; i < tensorSize; ++i) {
      devStrainTpdt = strain[i] - diag[i]*meanStrainTpdt;
      devStrainT = stateVars[s_totalStrain+i] - diag[i]*meanStrainT;
      stateVarsE[s_viscousStrain+i] = 
	expFac * stateVars[s_viscousStrain+i] + 
	dq * (devStrainTpdt - devStrainT);
    } //for

    material._updateStateVars(&stateVars[0], stateVars.size(), 
			      &properties[0], properties.size(),
			      &strain[0], strain.size(),
			      &initialStress[0], initialStress.size(),
			      &initialStrain[0], initialStrain.size());
    
    const double tolerance = 1.0e-06;
    for (int i=0; i < numVarsQuadPt; ++i)
      CPPUNIT_ASSERT_DOUBLES_EQUAL(stateVarsE[i], stateVars[i], tolerance);
  } // for




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
    initialState[i] = 0.1*i;
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
  
  material._updateStateVars(&parameters[0], numParamsQuadPt, 
			&totalStrainTpdt[0], tensorSize,
			&initialState[0], initialStateSize);

  const double tolerance = 1.0e-06;
  for (int i=0; i < numParamsQuadPt; ++i)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(parametersE[i], parameters[i], tolerance);
#endif
} // testUpdateStateVarsTimeDep

// ----------------------------------------------------------------------
// Test _stableTimeStepImplicit()
void
pylith::materials::TestGenMaxwellIsotropic3D::test_stableTimeStepImplicit(void)
{ // test_stableTimeStepImplicit
#if 0
  CPPUNIT_ASSERT(0 != _matElastic);

  delete _dataElastic; _dataElastic = new GenMaxwellIsotropic3DTimeDepData();

  TestElasticMaterial::test_stableTimeStepImplicit();
#endif
} // test_stableTimeStepImplicit


// End of file 
