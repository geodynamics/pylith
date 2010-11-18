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
// Copyright (c) 2010 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "TestGenMaxwellQpQsIsotropic3D.hh" // Implementation of class methods

#include "data/GenMaxwellQpQsIsotropic3DElasticData.hh" // USES GenMaxwellQpQsIsotropic3DElasticData
#include "data/GenMaxwellQpQsIsotropic3DTimeDepData.hh" // USES GenMaxwellQpQsIsotropic3DTimeDepData

#include "pylith/materials/GenMaxwellQpQsIsotropic3D.hh" // USES GenMaxwellQpQsIsotropic3D

#include <cstring> // USES memcpy()

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::materials::TestGenMaxwellQpQsIsotropic3D );

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::materials::TestGenMaxwellQpQsIsotropic3D::setUp(void)
{ // setUp
  _material = new GenMaxwellQpQsIsotropic3D();
  _matElastic = new GenMaxwellQpQsIsotropic3D();

  _data = new GenMaxwellQpQsIsotropic3DElasticData();
  _dataElastic = new GenMaxwellQpQsIsotropic3DElasticData();
  setupNormalizer();
} // setUp

// ----------------------------------------------------------------------
// Test timeStep()
void
pylith::materials::TestGenMaxwellQpQsIsotropic3D::testTimeStep(void)
{ // testTimeStep
  GenMaxwellQpQsIsotropic3D material;

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
pylith::materials::TestGenMaxwellQpQsIsotropic3D::testUseElasticBehavior(void)
{ // testUseElasticBehavior
  GenMaxwellQpQsIsotropic3D material;

  material.useElasticBehavior(true);
  // Some compilers/operating systems (cygwin) don't allow comparing
  // pointers. Use first test to determine if we can compare pointers.
  if (&pylith::materials::GenMaxwellQpQsIsotropic3D::_calcStressElastic == 
      material._calcStressFn) {
    CPPUNIT_ASSERT(&pylith::materials::GenMaxwellQpQsIsotropic3D::_calcStressElastic ==
		   material._calcStressFn);
    CPPUNIT_ASSERT(&pylith::materials::GenMaxwellQpQsIsotropic3D::_calcElasticConstsElastic ==
		   material._calcElasticConstsFn);
    CPPUNIT_ASSERT(&pylith::materials::GenMaxwellQpQsIsotropic3D::_updateStateVarsElastic ==
		   material._updateStateVarsFn);

    material.useElasticBehavior(false);
    CPPUNIT_ASSERT(&pylith::materials::GenMaxwellQpQsIsotropic3D::_calcStressViscoelastic ==
		   material._calcStressFn);
    CPPUNIT_ASSERT(&pylith::materials::GenMaxwellQpQsIsotropic3D::_calcElasticConstsViscoelastic == 
		   material._calcElasticConstsFn);
  CPPUNIT_ASSERT(&pylith::materials::GenMaxwellQpQsIsotropic3D::_updateStateVarsViscoelastic ==
		 material._updateStateVarsFn);
  } // if
} // testUseElasticBehavior

// ----------------------------------------------------------------------
// Test usesHasStateVars()
void
pylith::materials::TestGenMaxwellQpQsIsotropic3D::testHasStateVars(void)
{ // testHasStateVars
  GenMaxwellQpQsIsotropic3D material;
  CPPUNIT_ASSERT_EQUAL(true, material.hasStateVars());
} // testHasStateVars

// ----------------------------------------------------------------------
// Test calcStressElastic()
void
pylith::materials::TestGenMaxwellQpQsIsotropic3D::test_calcStressElastic(void)
{ // testCalcStressElastic
  CPPUNIT_ASSERT(0 != _matElastic);
  _matElastic->useElasticBehavior(true);

  test_calcStress();
} // testCalcStressElastic

// ----------------------------------------------------------------------
// Test calcElasticConstsElastic()
void
pylith::materials::TestGenMaxwellQpQsIsotropic3D::test_calcElasticConstsElastic(void)
{ // testElasticConstsElastic
  CPPUNIT_ASSERT(0 != _matElastic);
  _matElastic->useElasticBehavior(true);

  test_calcElasticConsts();
} // testElasticConstsElastic

// ----------------------------------------------------------------------
// Test updateStateVarsElastic()
void
pylith::materials::TestGenMaxwellQpQsIsotropic3D::test_updateStateVarsElastic(void)
{ // testUpdateStateVarsElastic
  // :TODO: Use TestElasticMaterial::test_updateStateVars
  // instead. This requires moving the calculation of the expected
  // state vars below to the Python code (where it belongs) and
  // setting the stateVarsUpdate attribute in the Python object.

  GenMaxwellQpQsIsotropic3D material;
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

    const double diag[] = { 1.0, 1.0, 1.0, 0.0, 0.0, 0.0 };
    const int numMaxwellModels = 3;

    const double meanStrain = (strain[0] + strain[1] + strain[2]) / 3.0;
    
    // Compute expected state variables
    double_array stateVarsE(numVarsQuadPt);
    const int s_totalStrain = 0;
    const int s_viscousStrain = s_totalStrain + tensorSize;
    const int s_viscousStrainBulk = s_viscousStrain + tensorSize*numMaxwellModels;

    // State variable 'total_strain' should match 'strain'
    for (int i=0; i < tensorSize; ++i) 
      stateVarsE[s_totalStrain+i] = strain[i];
    
    // State variable 'viscous_strain'
    for (int imodel=0; imodel < numMaxwellModels; ++imodel)
      for (int i=0; i < tensorSize; ++i)
	stateVarsE[s_viscousStrain+imodel*tensorSize+i] =
	  strain[i] - diag[i]*meanStrain;

    for (int imodel=0; imodel < numMaxwellModels; ++imodel)
      stateVarsE[s_viscousStrainBulk+imodel] = meanStrain;

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
pylith::materials::TestGenMaxwellQpQsIsotropic3D::test_calcStressTimeDep(void)
{ // testCalcStressTimeDep
  CPPUNIT_ASSERT(0 != _matElastic);
  _matElastic->useElasticBehavior(false);

  delete _dataElastic; _dataElastic = new GenMaxwellQpQsIsotropic3DTimeDepData();

  double dt = 2.0e+5;
  _matElastic->timeStep(dt);
  test_calcStress();
} // testCalcStressTimeDep

// ----------------------------------------------------------------------
// Test calcElasticConstsTimeDep()
void
pylith::materials::TestGenMaxwellQpQsIsotropic3D::test_calcElasticConstsTimeDep(void)
{ // testElasticConstsTimeDep
  CPPUNIT_ASSERT(0 != _matElastic);
  _matElastic->useElasticBehavior(false);

  delete _dataElastic; _dataElastic = new GenMaxwellQpQsIsotropic3DTimeDepData();

  double dt = 2.0e+5;
  _matElastic->timeStep(dt);
  test_calcElasticConsts();
} // testElasticConstsTimeDep

// ----------------------------------------------------------------------
// Test updateStateVarsTimeDep()
void
pylith::materials::TestGenMaxwellQpQsIsotropic3D::test_updateStateVarsTimeDep(void)
{ // testUpdateStateVarsTimeDep
  // :TODO: Use TestElasticMaterial::test_updateStateVars
  // instead. This requires moving the calculation of the expected
  // state vars below to the Python code (where it belongs) and
  // setting the stateVarsUpdate attribute in the Python object.

  GenMaxwellQpQsIsotropic3D material;
  material.useElasticBehavior(false);

  delete _dataElastic; _dataElastic = new GenMaxwellQpQsIsotropic3DTimeDepData();

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

    // Compute expected state variables
    double_array stateVarsE(numVarsQuadPt);
    const int numMaxwellModels = 3;
    const int s_totalStrain = 0;
    const int s_viscousStrain = s_totalStrain + tensorSize;
    const int s_viscousStrainBulk = s_viscousStrain + tensorSize*numMaxwellModels;
    const int p_shearRatio = 3;
    const int p_bulkRatio = p_shearRatio + numMaxwellModels;
    const int p_maxwellTimeShear = p_bulkRatio + numMaxwellModels;
    const int p_maxwellTimeBulk = p_maxwellTimeShear + numMaxwellModels;

    // State variable 'total_strain' should match 'strain'
    for (int i=0; i < tensorSize; ++i) 
      stateVarsE[s_totalStrain+i] = strain[i];
    
    // State variable 'viscous_strain'
    double_array maxwellTime(numMaxwellModels);
    double_array maxwellTimeBulk(numMaxwellModels);
    //    double_array shearRatio(numMaxwellModels);
    double_array dq(numMaxwellModels);
    double_array dqBulk(numMaxwellModels);
    for (int i=0; i < numMaxwellModels; ++i) {
      //shearRatio[i] = properties[p_shearRatio+i];
      maxwellTime[i] = properties[p_maxwellTimeShear+i];
      dq[i] = maxwellTime[i] * (1.0 - exp(-dt/maxwellTime[i]))/dt;
    } // for
    for (int i=0; i < numMaxwellModels; ++i) {
      //shearRatio[i] = properties[p_bulkRatio+i];
      maxwellTimeBulk[i] = properties[p_maxwellTimeBulk+i];
      dqBulk[i] = maxwellTimeBulk[i] * (1.0 - exp(-dt/maxwellTimeBulk[i]))/dt;
    } // for
    double_array strainT(tensorSize);
    for (int i=0; i < tensorSize; ++i)
      strainT[i] = stateVars[s_totalStrain+i];
    const double meanStrainT = 
      (stateVars[s_totalStrain+0] + 
       stateVars[s_totalStrain+1] + 
       stateVars[s_totalStrain+2]) / 3.0;
    const double meanStrainTpdt = (strain[0] + strain[1] + strain[2]) / 3.0;

    double devStrainTpdt = 0.0;
    double devStrainT = 0.0;
    double deltaStrain = 0.0;
    double visStrain = 0.0;
    const double diag[] = { 1.0, 1.0, 1.0, 0.0, 0.0, 0.0 };
    
    for (int iComp=0; iComp < tensorSize; ++iComp) {
      devStrainTpdt = strain[iComp] - diag[iComp]*meanStrainTpdt;
      devStrainT = strainT[iComp] - diag[iComp]*meanStrainT;
      deltaStrain = devStrainTpdt - devStrainT;
      for (int imodel=0; imodel < numMaxwellModels; ++imodel) {
	stateVarsE[s_viscousStrain+imodel*tensorSize+iComp] =
	  exp(-dt/maxwellTime[imodel]) *
	  stateVars[s_viscousStrain+imodel*tensorSize+iComp] + dq[imodel] * deltaStrain;
      } // for
    } // for

    double volStrainTpdt = meanStrainTpdt;
    double volStrainT = meanStrainT;
    deltaStrain = volStrainTpdt - volStrainT;

    for (int imodel=0; imodel < numMaxwellModels; ++imodel) {
      stateVarsE[s_viscousStrainBulk+imodel] = 
	  exp(-dt/maxwellTimeBulk[imodel]) *
	  stateVars[s_viscousStrainBulk+imodel] + dqBulk[imodel] * deltaStrain;
      } // for

    material._updateStateVars(&stateVars[0], stateVars.size(), 
			      &properties[0], properties.size(),
			      &strain[0], strain.size(),
			      &initialStress[0], initialStress.size(),
			      &initialStrain[0], initialStrain.size());

    const double tolerance = 1.0e-06;
    for (int i=0; i < numVarsQuadPt; ++i)
      if (stateVarsE[i] > tolerance)
	CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, stateVars[i]/stateVarsE[i], tolerance);
      else
	CPPUNIT_ASSERT_DOUBLES_EQUAL(stateVarsE[i], stateVars[i], tolerance);
  } // for
} // testUpdateStateVarsTimeDep

// ----------------------------------------------------------------------
// Test _stableTimeStepImplicit()
void
pylith::materials::TestGenMaxwellQpQsIsotropic3D::test_stableTimeStepImplicit(void)
{ // test_stableTimeStepImplicit
  CPPUNIT_ASSERT(0 != _matElastic);

  delete _dataElastic; _dataElastic = new GenMaxwellQpQsIsotropic3DTimeDepData();

  TestElasticMaterial::test_stableTimeStepImplicit();
} // test_stableTimeStepImplicit


// End of file 
