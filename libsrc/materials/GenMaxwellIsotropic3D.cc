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

#include "GenMaxwellIsotropic3D.hh" // implementation of object methods

#include "ViscoelasticMaxwell.hh" // USES computeVisStrain

#include "pylith/utils/array.hh" // USES double_array

#include "petsc.h" // USES PetscLogFlopsNoCheck

#include <assert.h> // USES assert()
#include <sstream> // USES std::ostringstream
#include <stdexcept> // USES std::runtime_error

// ----------------------------------------------------------------------
namespace pylith {
  namespace materials {
    namespace _GenMaxwellIsotropic3D{

      /// Number of entries in stress/strain tensors.
      const int tensorSize = 6;

      /// Number of entries in derivative of elasticity matrix.
      const int numElasticConsts = 21;

      /// Number of Maxwell models in parallel.
      const int numMaxwellModels = 3;

      /// Values expected in spatial database
      const int numDBValues = 3 + 2*numMaxwellModels;
      /// NOTE:  I haven't generalized the database names to the number
      /// of Maxwell models like I have for everything else.
      const char* namesDBValues[] = {"density", "vs", "vp" ,
				     "shear_ratio_1",
				     "shear_ratio_2",
				     "shear_ratio_3",
				     "viscosity_1",
				     "viscosity_2",
				     "viscosity_3"};

      /// Indices (order) of database values
      const int didDensity = 0;
      const int didVs = 1;
      const int didVp = 2;
      const int didShearRatio1 = 3;
      const int didShearRatio2 = 4;
      const int didShearRatio3 = 5;
      const int didViscosity1 = 6;
      const int didViscosity2 = 7;
      const int didViscosity3 = 8;

      /// Parameters
      const int numParameters = 7;
      const int numParamValues[] = { 1, 
				     1, 
				     1,
				     numMaxwellModels,
				     numMaxwellModels,
				     tensorSize,
				     tensorSize * numMaxwellModels};

      /// Indices (order) of parameters
      const int pidDensity = 0;
      const int pidMuTot = pidDensity + 1;
      const int pidLambdaTot = pidMuTot + 1;
      const int pidShearRatio = pidLambdaTot + 1;
      const int pidMaxwellTime = pidShearRatio + numMaxwellModels;
      const int pidStrainT = pidMaxwellTime + numMaxwellModels;
      const int pidVisStrain = pidStrainT + tensorSize;
    } // _GenMaxwellIsotropic3D
  } // materials
} // pylith

// ----------------------------------------------------------------------
// Default constructor.
pylith::materials::GenMaxwellIsotropic3D::GenMaxwellIsotropic3D(void) :
  ElasticMaterial(_GenMaxwellIsotropic3D::numParamValues,
		  _GenMaxwellIsotropic3D::numParameters),
  _calcElasticConstsFn(&pylith::materials::GenMaxwellIsotropic3D::_calcElasticConstsElastic),
  _calcStressFn(&pylith::materials::GenMaxwellIsotropic3D::_calcStressElastic),
  _updateStateFn(&pylith::materials::GenMaxwellIsotropic3D::_updateStateElastic)
{ // constructor
  _dimension = 3;
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::materials::GenMaxwellIsotropic3D::~GenMaxwellIsotropic3D(void)
{ // destructor
} // destructor

// ----------------------------------------------------------------------
// Get names of values expected to be in database of parameters for
const char**
pylith::materials::GenMaxwellIsotropic3D::_dbValues(void) const
{ // _dbValues
  return _GenMaxwellIsotropic3D::namesDBValues;
} // _dbValues

// ----------------------------------------------------------------------
// Get number of values expected to be in database of parameters for
int
pylith::materials::GenMaxwellIsotropic3D::_numDBValues(void) const
{ // _numDBValues
  return _GenMaxwellIsotropic3D::numDBValues;
} // _numDBValues

// ----------------------------------------------------------------------
// Compute parameters from values in spatial database.
void
pylith::materials::GenMaxwellIsotropic3D::_dbToParameters(
					    double* const paramVals,
					    const int numParams,
					    const double_array& dbValues) const
{ // _dbToParameters
  assert(0 != paramVals);
  assert(_numParamsQuadPt == numParams);
  const int numDBValues = dbValues.size();
  assert(_GenMaxwellIsotropic3D::numDBValues == numDBValues);

  const double density = dbValues[_GenMaxwellIsotropic3D::didDensity];
  const double vs = dbValues[_GenMaxwellIsotropic3D::didVs];
  const double vp = dbValues[_GenMaxwellIsotropic3D::didVp];
  const int numMaxwellModels = _GenMaxwellIsotropic3D::numMaxwellModels;
 
  const double mu = density * vs*vs;
  const double lambda = density * vp*vp - 2.0*mu;
  if (lambda < 0.0) {
    std::ostringstream msg;
    msg << "Attempted to set Lame's constant lambda to negative value.\n"
	<< "density: " << density << "\n"
	<< "vp: " << vp << "\n"
	<< "vs: " << vs << "\n";
    throw std::runtime_error(msg.str());
  } // if

  assert(mu > 0);

  paramVals[_GenMaxwellIsotropic3D::pidDensity] = density;
  paramVals[_GenMaxwellIsotropic3D::pidMuTot] = mu;
  paramVals[_GenMaxwellIsotropic3D::pidLambdaTot] = lambda;

  // Loop over number of Maxwell models.
  for (int iTerm =0;
       iTerm < numMaxwellModels; ++iTerm) {
    double muRatio = dbValues[_GenMaxwellIsotropic3D::didShearRatio1 + iTerm];
    double viscosity = dbValues[_GenMaxwellIsotropic3D::didViscosity1 + iTerm];
    double mu = muRatio*mu;
    double maxwellTime = 0.0;
    if (mu > 0.0) maxwellTime = viscosity / mu;
    paramVals[_GenMaxwellIsotropic3D::pidShearRatio + iTerm] = muRatio;
    paramVals[_GenMaxwellIsotropic3D::pidMaxwellTime + iTerm] = maxwellTime;
  } // for

  PetscLogFlopsNoCheck(6+2*numMaxwellModels);
} // _dbToParameters

// ----------------------------------------------------------------------
// Get number of entries in stress tensor.
int
pylith::materials::GenMaxwellIsotropic3D::_tensorSize(void) const
{ // _tensorSize
  return _GenMaxwellIsotropic3D::tensorSize;
} // _tensorSize

// ----------------------------------------------------------------------
// Get number of entries in elasticity matrix for material.
int
pylith::materials::GenMaxwellIsotropic3D::_numElasticConsts(void) const
{ // numElasticConsts
  return _GenMaxwellIsotropic3D::numElasticConsts;
} // numElasticConsts

// ----------------------------------------------------------------------
// Compute density at location from parameters.
void
pylith::materials::GenMaxwellIsotropic3D::_calcDensity(double* const density,
						       const double* parameters,
						       const int numParams)
{ // _calcDensity
  assert(0 != density);
  assert(0 != parameters);
  assert(_numParamsQuadPt == numParams);
  
  density[0] = parameters[_GenMaxwellIsotropic3D::pidDensity];
} // _calcDensity

// ----------------------------------------------------------------------
// Compute stress tensor at location from parameters as an elastic
// material.
void
pylith::materials::GenMaxwellIsotropic3D::_calcStressElastic(
							     double* const stress,
							     const int stressSize,
							     const double* parameters,
							     const int numParams,
							     const double* totalStrain,
							     const int strainSize)
{ // _calcStressElastic
  assert(0 != stress);
  assert(_GenMaxwellIsotropic3D::tensorSize == stressSize);
  assert(0 != parameters);
  assert(_numParamsQuadPt == numParams);
  assert(0 != totalStrain);
  assert(_GenMaxwellIsotropic3D::tensorSize == strainSize);

  const double density = parameters[_GenMaxwellIsotropic3D::pidDensity];
  const double mu = parameters[_GenMaxwellIsotropic3D::pidMuTot];
  const double lambda = parameters[_GenMaxwellIsotropic3D::pidLambdaTot];

  const double mu2 = 2.0 * mu;

  const double e11 = totalStrain[0];
  const double e22 = totalStrain[1];
  const double e33 = totalStrain[2];
  const double e12 = totalStrain[3];
  const double e23 = totalStrain[4];
  const double e13 = totalStrain[5];
  
  const double traceStrainTpdt = e11 + e22 + e33;
  const double s123 = lambda * traceStrainTpdt;

  stress[0] = s123 + mu2*e11;
  stress[1] = s123 + mu2*e22;
  stress[2] = s123 + mu2*e33;
  stress[3] = mu2 * e12;
  stress[4] = mu2 * e23;
  stress[5] = mu2 * e13;
  // std::cout << " _calcStressElastic: " << std::endl;
  // std::cout << " totalStrain: " << std::endl;
  // for (int iComp=0; iComp < _GenMaxwellIsotropic3D::tensorSize; ++iComp)
    // std::cout << "  " << totalStrain[iComp];
  // std::cout << std::endl << " stress: " << std::endl;
  // for (int iComp=0; iComp < _GenMaxwellIsotropic3D::tensorSize; ++iComp)
    // std::cout << "  " << (*stress)[iComp];
  // std::cout << std::endl;

  PetscLogFlopsNoCheck(13);
} // _calcStressElastic

// ----------------------------------------------------------------------
// Compute stress tensor at location from parameters as a viscoelastic
// material.
void
pylith::materials::GenMaxwellIsotropic3D::_calcStressViscoelastic(
								  double* const stress,
								  const int stressSize,
								  const double* parameters,
								  const int numParams,
								  const double* totalStrain,
								  const int strainSize)
{ // _calcStressElastic
  assert(0 != stress);
  assert(_GenMaxwellIsotropic3D::tensorSize == stressSize);
  assert(0 != parameters);
  assert(_numParamsQuadPt == numParams);
  assert(0 != totalStrain);
  assert(_GenMaxwellIsotropic3D::tensorSize == strainSize);

  const double density = parameters[_GenMaxwellIsotropic3D::pidDensity];
  const double mu = parameters[_GenMaxwellIsotropic3D::pidMuTot];
  const double lambda = parameters[_GenMaxwellIsotropic3D::pidLambdaTot];
  const double muRatio[] = {
    parameters[_GenMaxwellIsotropic3D::pidShearRatio],
    parameters[_GenMaxwellIsotropic3D::pidShearRatio+1],
    parameters[_GenMaxwellIsotropic3D::pidShearRatio+2]};
  const double maxwellTime[] = {
    parameters[_GenMaxwellIsotropic3D::pidMaxwellTime],
    parameters[_GenMaxwellIsotropic3D::pidMaxwellTime+1],
    parameters[_GenMaxwellIsotropic3D::pidMaxwellTime+2]};
  const int numMaxwellModels = _GenMaxwellIsotropic3D::numMaxwellModels;

  const double mu2 = 2.0 * mu;
  const double bulkmodulus = lambda + mu2/3.0;

  const double e11 = totalStrain[0];
  const double e22 = totalStrain[1];
  const double e33 = totalStrain[2];
  const double e12 = totalStrain[3];
  const double e23 = totalStrain[4];
  const double e13 = totalStrain[5];
  
  const double traceStrainTpdt = e11 + e22 + e33;
  const double meanStrainTpdt = traceStrainTpdt/3.0;
  const double s123 = lambda * traceStrainTpdt;

  const double diag[] = { 1.0, 1.0, 1.0, 0.0, 0.0, 0.0 };
  const double meanStressTpdt = bulkmodulus * traceStrainTpdt;

  const double meanStrainT = 
    (parameters[_GenMaxwellIsotropic3D::pidStrainT+0] +
     parameters[_GenMaxwellIsotropic3D::pidStrainT+1] +
     parameters[_GenMaxwellIsotropic3D::pidStrainT+2])/3.0;
  
  double visFrac = 0.0;
  for (int iTerm = 0; iTerm < numMaxwellModels; ++iTerm) 
    visFrac += muRatio[iTerm];

  if (visFrac > 1.0) {
    std::ostringstream msg;
    msg << "Shear modulus ratios sum to a value greater than 1.0 for\n"
	<< "Generalized Maxwell model.\n"
	<< "Ratio 1: " << muRatio[0] << "\n"
	<< "Ratio 2: " << muRatio[1] << "\n"
	<< "Ratio 3: " << muRatio[2] << "\n"
	<< "Total:   " << visFrac << "\n";
    throw std::runtime_error(msg.str());
  } // if
  double elasFrac = 1.0 - visFrac;

  PetscLogFlopsNoCheck(12+numMaxwellModels);

  // Compute Prony series terms
  double dq[] = { 0.0, 0.0, 0.0};
  for (int iTerm=0; iTerm < numMaxwellModels; ++iTerm) {
    if (muRatio[iTerm] != 0.0)
      dq[iTerm] =
	ViscoelasticMaxwell::computeVisStrain(_dt, maxwellTime[iTerm]);
  } // for

  // Compute new stresses
  double devStrainTpdt = 0.0;
  double devStrainT = 0.0;
  double deltaStrain = 0.0;
  double devStressTpdt = 0.0;
  double visStrain = 0.0;
  // std::cout << " _calcStressViscoelastic: " << std::endl;
  // std::cout << " stress  totalStrain  visStrain: " << std::endl;
  for (int iComp=0; iComp < _GenMaxwellIsotropic3D::tensorSize; ++iComp) {
    devStrainTpdt = totalStrain[iComp] - diag[iComp]*meanStrainTpdt;
    devStrainT = parameters[_GenMaxwellIsotropic3D::pidStrainT+iComp] -
      diag[iComp]*meanStrainT;
    deltaStrain = devStrainTpdt - devStrainT;
    devStressTpdt = elasFrac*devStrainTpdt;
    for (int iTerm=0; iTerm < numMaxwellModels; ++iTerm) {
      if (muRatio[iTerm] != 0.0) {
	visStrain = exp(-_dt/maxwellTime[iTerm])*
	  parameters[_GenMaxwellIsotropic3D::pidVisStrain+
		     iComp+iTerm*_GenMaxwellIsotropic3D::tensorSize]
	  + dq[iTerm]*deltaStrain;
	devStressTpdt += muRatio[iTerm] * visStrain;
      } // if
    } // for

    devStressTpdt = mu2*devStressTpdt;
    // Later I will want to put in initial stresses.
    stress[iComp] =diag[iComp]*meanStressTpdt+devStressTpdt;

    // Temporary to get stresses and strains.
    // std::cout << "  " << stress[iComp] << "  " << totalStrain[iComp] << "  " << visStrain << std:: endl;
  } // for

  PetscLogFlopsNoCheck((8 + 8*numMaxwellModels) *
		       _GenMaxwellIsotropic3D::tensorSize);
} // _calcStress

// ----------------------------------------------------------------------
// Compute derivative of elasticity matrix at location from parameters.
void
pylith::materials::GenMaxwellIsotropic3D::_calcElasticConstsElastic(
						  double* const elasticConsts,
						  const int numElasticConsts,
						  const double* parameters,
						  const int numParams,
						  const double* totalStrain,
						  const int strainSize)
{ // _calcElasticConstsElastic
  assert(0 != elasticConsts);
  assert(_GenMaxwellIsotropic3D::numElasticConsts == numElasticConsts);
  assert(0 != parameters);
  assert(_numParamsQuadPt == numParams);
  assert(0 != totalStrain);
  assert(_GenMaxwellIsotropic3D::tensorSize == strainSize);
 
  const double density = parameters[_GenMaxwellIsotropic3D::pidDensity];
  const double mu = parameters[_GenMaxwellIsotropic3D::pidMuTot];
  const double lambda = parameters[_GenMaxwellIsotropic3D::pidLambdaTot];

  const double mu2 = 2.0 * mu;
  const double lambda2mu = lambda + mu2;
  const double bulkmodulus = lambda + mu2 / 3.0;

  elasticConsts[ 0] = lambda2mu; // C1111
  elasticConsts[ 1] = lambda; // C1122
  elasticConsts[ 2] = lambda; // C1133
  elasticConsts[ 3] = 0; // C1112
  elasticConsts[ 4] = 0; // C1123
  elasticConsts[ 5] = 0; // C1113
  elasticConsts[ 6] = lambda2mu; // C2222
  elasticConsts[ 7] = lambda; // C2233
  elasticConsts[ 8] = 0; // C2212
  elasticConsts[ 9] = 0; // C2223
  elasticConsts[10] = 0; // C2213
  elasticConsts[11] = lambda2mu; // C3333
  elasticConsts[12] = 0; // C3312
  elasticConsts[13] = 0; // C3323
  elasticConsts[14] = 0; // C3313
  elasticConsts[15] = mu2; // C1212
  elasticConsts[16] = 0; // C1223
  elasticConsts[17] = 0; // C1213
  elasticConsts[18] = mu2; // C2323
  elasticConsts[19] = 0; // C2313
  elasticConsts[20] = mu2; // C1313

  PetscLogFlopsNoCheck(4);
} // _calcElasticConstsElastic

// ----------------------------------------------------------------------
// Compute derivative of elasticity matrix at location from parameters
// as an elastic material.
void
pylith::materials::GenMaxwellIsotropic3D::_calcElasticConstsViscoelastic(
						  double* const elasticConsts,
						  const int numElasticConsts,
						  const double* parameters,
						  const int numParams,
						  const double* totalStrain,
						  const int strainSize)
{ // _calcElasticConstsViscoelastic
  assert(0 != elasticConsts);
  assert(_GenMaxwellIsotropic3D::numElasticConsts == numElasticConsts);
  assert(0 != parameters);
  assert(_numParamsQuadPt == numParams);
  assert(0 != totalStrain);
  assert(_GenMaxwellIsotropic3D::tensorSize == strainSize);
 
  const double density = parameters[_GenMaxwellIsotropic3D::pidDensity];
  const double mu = parameters[_GenMaxwellIsotropic3D::pidMuTot];
  const double lambda = parameters[_GenMaxwellIsotropic3D::pidLambdaTot];
  const int numMaxwellModels = _GenMaxwellIsotropic3D::numMaxwellModels;

  const double mu2 = 2.0 * mu;
  const double bulkmodulus = lambda + mu2/3.0;

  // Compute viscous contribution.
  double visFac = 0.0;
  double shearRatio = 0.0;
  double maxwellTime = 0.0;
  for (int iTerm = 0; iTerm < numMaxwellModels; ++iTerm) {
    shearRatio = parameters[_GenMaxwellIsotropic3D::pidShearRatio + iTerm];
    if (shearRatio != 0.0) {
      maxwellTime = parameters[_GenMaxwellIsotropic3D::pidMaxwellTime + iTerm];
      visFac +=
	shearRatio*ViscoelasticMaxwell::computeVisStrain(_dt, maxwellTime);
    } // if
  } // for
  visFac = mu*visFac/3.0;

  elasticConsts[ 0] = bulkmodulus + 4.0*visFac; // C1111
  elasticConsts[ 1] = bulkmodulus - 2.0*visFac; // C1122
  elasticConsts[ 2] = elasticConsts[1]; // C1133
  elasticConsts[ 3] = 0; // C1112
  elasticConsts[ 4] = 0; // C1123
  elasticConsts[ 5] = 0; // C1113
  elasticConsts[ 6] = elasticConsts[0]; // C2222
  elasticConsts[ 7] = elasticConsts[1]; // C2233
  elasticConsts[ 8] = 0; // C2212
  elasticConsts[ 9] = 0; // C2223
  elasticConsts[10] = 0; // C2213
  elasticConsts[11] = elasticConsts[0]; // C3333
  elasticConsts[12] = 0; // C3312
  elasticConsts[13] = 0; // C3323
  elasticConsts[14] = 0; // C3313
  elasticConsts[15] = 6.0 * visFac; // C1212
  elasticConsts[16] = 0; // C1223
  elasticConsts[17] = 0; // C1213
  elasticConsts[18] = elasticConsts[15]; // C2323
  elasticConsts[19] = 0; // C2313
  elasticConsts[20] = elasticConsts[15]; // C1313

  PetscLogFlopsNoCheck(8 + 2*numMaxwellModels);
} // _calcElasticConstsViscoelastic

// ----------------------------------------------------------------------
// Update state variables.
void
pylith::materials::GenMaxwellIsotropic3D::_updateStateElastic(
						 double* const parameters,
						 const int numParams,
						 const double* totalStrain,
						 const int strainSize)
{ // _updateStateElastic
  assert(0 != parameters);
  assert(_numParamsQuadPt == numParams);
  assert(0 != totalStrain);
  assert(_GenMaxwellIsotropic3D::tensorSize == strainSize);

  const int numMaxwellModels = _GenMaxwellIsotropic3D::numMaxwellModels;

  const double e11 = totalStrain[0];
  const double e22 = totalStrain[1];
  const double e33 = totalStrain[2];

  const double traceStrainTpdt = e11 + e22 + e33;
  const double meanStrainTpdt = traceStrainTpdt/3.0;

  const double diag[] = { 1.0, 1.0, 1.0, 0.0, 0.0, 0.0 };

  // Temporary to get stresses.
  // double_array stress(6);
  // _calcStressElastic(&stress, (*parameters), totalStrain);

  // Initialize all viscous strains to deviatoric elastic strains.
  double devStrain = 0.0;
  for (int iComp=0; iComp < _GenMaxwellIsotropic3D::tensorSize; ++iComp) {
    parameters[_GenMaxwellIsotropic3D::pidStrainT+iComp] = totalStrain[iComp];
    devStrain = totalStrain[iComp] - diag[iComp]*meanStrainTpdt;
    for (int iTerm = 0; iTerm < numMaxwellModels; ++iTerm) {
      parameters[_GenMaxwellIsotropic3D::pidVisStrain+
		 iComp+iTerm*_GenMaxwellIsotropic3D::tensorSize] =
	devStrain;
    } // for
  } // for
  PetscLogFlopsNoCheck(3+2*_GenMaxwellIsotropic3D::tensorSize);
//   std::cout << std::endl;
//   std::cout << " updateStateElastic: "<< std::endl;
//   std::cout << " StrainT  VisStrain  Stress: " << std::endl;
//   for (int iComp=0; iComp < _GenMaxwellIsotropic3D::tensorSize; ++iComp)
//     std::cout << "  " << parameters[_GenMaxwellIsotropic3D::pidStrainT+iComp]
// 	    << "   " << parameters[_GenMaxwellIsotropic3D::pidVisStrain+iComp]
// 	    << "   " << stress[iComp]
// 	    << std::endl;

  _needNewJacobian = true;
} // _updateStateElastic

// ----------------------------------------------------------------------
// Update state variables.
void
pylith::materials::GenMaxwellIsotropic3D::_updateStateViscoelastic(
						 double* const parameters,
						 const int numParams,
						 const double* totalStrain,
						 const int strainSize)
{ // _updateStateViscoelastic
  assert(0 != parameters);
  assert(_numParamsQuadPt == numParams);
  assert(0 != totalStrain);
  assert(_GenMaxwellIsotropic3D::tensorSize == strainSize);

  const int numMaxwellModels = _GenMaxwellIsotropic3D::numMaxwellModels;
  const int tensorSize = _GenMaxwellIsotropic3D::tensorSize;

  const double e11 = totalStrain[0];
  const double e22 = totalStrain[1];
  const double e33 = totalStrain[2];

  const double meanStrainTpdt = (e11 + e22 + e33) / 3.0;

  const double diag[] = { 1.0, 1.0, 1.0, 0.0, 0.0, 0.0 };

  // Temporary to get stresses.
  // double stress[6];
  // const int stressSize = 6;
  // _calcStressViscoelastic(stress, stressSize, 
			  // parameters, numParams,
			  // totalStrain, strainSize);

  const double meanStrainT = 
    (parameters[_GenMaxwellIsotropic3D::pidStrainT+0] +
     parameters[_GenMaxwellIsotropic3D::pidStrainT+1] +
     parameters[_GenMaxwellIsotropic3D::pidStrainT+2]) / 3.0;
  
  PetscLogFlopsNoCheck(6);

  // Compute Prony series terms.
  double dq[] = {0.0, 0.0, 0.0};
  for (int iTerm = 0; iTerm < numMaxwellModels; ++iTerm) {
    if(parameters[_GenMaxwellIsotropic3D::pidShearRatio + iTerm] != 0.0) {
      const double maxwellTime =
	parameters[_GenMaxwellIsotropic3D::pidMaxwellTime + iTerm];
      dq[iTerm] = ViscoelasticMaxwell::computeVisStrain(_dt, maxwellTime);
    } // if
  } // for
      
  double devStrainTpdt = 0.0;
  double devStrainT = 0.0;
  double deltaStrain = 0.0;
  double visStrain = 0.0;
  for (int iComp=0; iComp < tensorSize; ++iComp) {
    devStrainTpdt = totalStrain[iComp] - diag[iComp]*meanStrainTpdt;
    devStrainT = parameters[_GenMaxwellIsotropic3D::pidStrainT+iComp] -
      diag[iComp] * meanStrainT;
    deltaStrain = devStrainTpdt - devStrainT;
    parameters[_GenMaxwellIsotropic3D::pidStrainT+iComp] = totalStrain[iComp];
    for (int iTerm = 0; iTerm < numMaxwellModels; ++iTerm) {
      const double maxwellTime =
	parameters[_GenMaxwellIsotropic3D::pidMaxwellTime + iTerm];
      visStrain = 
	exp(-_dt/maxwellTime) * 
	parameters[_GenMaxwellIsotropic3D::pidVisStrain+iComp] +
	dq[iTerm] * deltaStrain;
      parameters[_GenMaxwellIsotropic3D::pidVisStrain +
		 iComp + iTerm * tensorSize] = visStrain;
    } // for
  } // for
  PetscLogFlopsNoCheck((5 + (6 * numMaxwellModels)) * tensorSize);

  _needNewJacobian = false;

//   std::cout << std::endl;
//   std::cout << " updateStateViscoelastic: "<< std::endl;
//   std::cout << " StrainT  VisStrain  Stress: " << std::endl;
//   for (int iComp=0; iComp < _GenMaxwellIsotropic3D::tensorSize; ++iComp)
//     std::cout << "  " << parameters[_GenMaxwellIsotropic3D::pidStrainT+iComp]
// 	    << "   " << parameters[_GenMaxwellIsotropic3D::pidVisStrain+iComp]
// 	    << "   " << stress[iComp]
// 	    << std::endl;
} // _updateStateViscoelastic


// End of file 
