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

#include <cstring> // USES memcpy()

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
  setupNormalizer();
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
  CPPUNIT_ASSERT_EQUAL(&pylith::materials::MaxwellIsotropic3D::_updateStateVarsElastic,
		       material._updateStateVarsFn);

  material.useElasticBehavior(false);
  CPPUNIT_ASSERT_EQUAL(&pylith::materials::MaxwellIsotropic3D::_calcStressViscoelastic,
		       material._calcStressFn);
  CPPUNIT_ASSERT_EQUAL(&pylith::materials::MaxwellIsotropic3D::_calcElasticConstsViscoelastic,
		       material._calcElasticConstsFn);
  CPPUNIT_ASSERT_EQUAL(&pylith::materials::MaxwellIsotropic3D::_updateStateVarsViscoelastic,
		       material._updateStateVarsFn);
} // testUseElasticBehavior

// ----------------------------------------------------------------------
// Test usesHasStateVars()
void
pylith::materials::TestMaxwellIsotropic3D::testHasStateVars(void)
{ // testHasStateVars
  MaxwellIsotropic3D material;
  CPPUNIT_ASSERT_EQUAL(true, material.hasStateVars());
} // testHasStateVars

// ----------------------------------------------------------------------
// Test _calcStressElastic()
void
pylith::materials::TestMaxwellIsotropic3D::test_calcStressElastic(void)
{ // test_calcStressElastic
  CPPUNIT_ASSERT(0 != _matElastic);
  _matElastic->useElasticBehavior(true);

  test_calcStress();
} // test_calcStressElastic

// ----------------------------------------------------------------------
// Test calcElasticConstsElastic()
void
pylith::materials::TestMaxwellIsotropic3D::test_calcElasticConstsElastic(void)
{ // test_calcElasticConstsElastic
  CPPUNIT_ASSERT(0 != _matElastic);
  _matElastic->useElasticBehavior(true);

  test_calcElasticConsts();
} // test_calcElasticConstsElastic

// ----------------------------------------------------------------------
// Test _updateStateVarsElastic()
void
pylith::materials::TestMaxwellIsotropic3D::test_updateStateVarsElastic(void)
{ // test_updateStateVarsElastic
  // :TODO: Use TestElasticMaterial::test_updateStateVars
  // instead. This requires moving the calculation of the expected
  // state vars below to the Python code (where it belongs) and
  // setting the stateVarsUpdate attribute in the Python object.

  MaxwellIsotropic3D material;
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
    for (int i=0; i < tensorSize; ++i)
      stateVarsE[s_viscousStrain+i] = strain[i] - diag[i]*meanStrain;
 
    material._updateStateVars(&stateVars[0], stateVars.size(), 
			      &properties[0], properties.size(),
			      &strain[0], strain.size(),
			      &initialStress[0], initialStress.size(),
			      &initialStrain[0], initialStrain.size());

    const double tolerance = 1.0e-06;
    for (int i=0; i < numVarsQuadPt; ++i)
      CPPUNIT_ASSERT_DOUBLES_EQUAL(stateVarsE[i], stateVars[i], tolerance);
  } // for
} // test_updateStateVarsElastic

// ----------------------------------------------------------------------
// Test _calcStressTimeDep()
void
pylith::materials::TestMaxwellIsotropic3D::test_calcStressTimeDep(void)
{ // test_calcStressTimeDep
  CPPUNIT_ASSERT(0 != _matElastic);
  _matElastic->useElasticBehavior(false);

  delete _dataElastic; _dataElastic = new MaxwellIsotropic3DTimeDepData();

  double dt = 2.0e+5;
  _matElastic->timeStep(dt);
  test_calcStress();
} // test_calcStressTimeDep

// ----------------------------------------------------------------------
// Test _calcElasticConstsTimeDep()
void
pylith::materials::TestMaxwellIsotropic3D::test_calcElasticConstsTimeDep(void)
{ // test_calcElasticConstsTimeDep
  CPPUNIT_ASSERT(0 != _matElastic);
  _matElastic->useElasticBehavior(false);

  delete _dataElastic; _dataElastic = new MaxwellIsotropic3DTimeDepData();

  double dt = 2.0e+5;
  _matElastic->timeStep(dt);
  test_calcElasticConsts();
} // test_calcElasticConstsTimeDep

// ----------------------------------------------------------------------
// Test _updateStateVarsTimeDep()
void
pylith::materials::TestMaxwellIsotropic3D::test_updateStateVarsTimeDep(void)
{ // test_updateStateVarsTimeDep
  // :TODO: Use TestElasticMaterial::test_updateStateVars
  // instead. This requires moving the calculation of the expected
  // state vars below to the Python code (where it belongs) and
  // setting the stateVarsUpdate attribute in the Python object.

  MaxwellIsotropic3D material;
  material.useElasticBehavior(false);
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
} // test_updateStateVarsTimeDep

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
