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

#include "MaxwellIsotropic3D.hh" // implementation of object methods

#include "ViscoelasticMaxwell.hh" // USES computeVisStrain

#include "pylith/utils/array.hh" // USES double_array

#include "petsc.h" // USES PetscLogFlopsNoCheck

#include <assert.h> // USES assert()
#include <sstream> // USES std::ostringstream
#include <stdexcept> // USES std::runtime_error

// ----------------------------------------------------------------------
namespace pylith {
  namespace materials {
    namespace _MaxwellIsotropic3D{

      /// Number of entries in stress/strain tensors.
      const int tensorSize = 6;

      /// Number of entries in derivative of elasticity matrix.
      const int numElasticConsts = 21;

      /// Values expected in spatial database
      const int numDBValues = 4;
      const char* namesDBValues[] = {"density", "vs", "vp" , "viscosity"};

      /// Indices (order) of database values
      const int didDensity = 0;
      const int didVs = 1;
      const int didVp = 2;
      const int didViscosity = 3;

      /// Parameters
      const int numParameters = 6;
      const int numParamValues[] = { 1, 1, 1, 1, 6, 6};

      /// Indices (order) of parameters
      const int pidDensity = 0;
      const int pidMu = pidDensity + 1;
      const int pidLambda = pidMu + 1;
      const int pidMaxwellTime = pidLambda + 1;
      const int pidStrainT = pidMaxwellTime + 1;
      const int pidVisStrain = pidStrainT + 6;
    } // _MaxwellIsotropic3D
  } // materials
} // pylith

// ----------------------------------------------------------------------
// Default constructor.
pylith::materials::MaxwellIsotropic3D::MaxwellIsotropic3D(void) :
  ElasticMaterial(_MaxwellIsotropic3D::numParamValues,
		  _MaxwellIsotropic3D::numParameters),
  _calcElasticConstsFn(&pylith::materials::MaxwellIsotropic3D::_calcElasticConstsElastic),
  _calcStressFn(&pylith::materials::MaxwellIsotropic3D::_calcStressElastic),
  _updateStateFn(&pylith::materials::MaxwellIsotropic3D::_updateStateElastic)
{ // constructor
  _dimension = 3;
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::materials::MaxwellIsotropic3D::~MaxwellIsotropic3D(void)
{ // destructor
} // destructor

// ----------------------------------------------------------------------
// Get names of values expected to be in database of parameters for
const char**
pylith::materials::MaxwellIsotropic3D::_dbValues(void) const
{ // _dbValues
  return _MaxwellIsotropic3D::namesDBValues;
} // _dbValues

// ----------------------------------------------------------------------
// Get number of values expected to be in database of parameters for
int
pylith::materials::MaxwellIsotropic3D::_numDBValues(void) const
{ // _numDBValues
  return _MaxwellIsotropic3D::numDBValues;
} // _numDBValues

// ----------------------------------------------------------------------
// Compute parameters from values in spatial database.
void
pylith::materials::MaxwellIsotropic3D::_dbToParameters(
					    double* const paramVals,
					    const int numParams,
					    const double_array& dbValues) const
{ // _dbToParameters
  assert(0 != paramVals);
  assert(_numParamsQuadPt == numParams);
  const int numDBValues = dbValues.size();
  assert(_MaxwellIsotropic3D::numDBValues == numDBValues);

  const double density = dbValues[_MaxwellIsotropic3D::didDensity];
  const double vs = dbValues[_MaxwellIsotropic3D::didVs];
  const double vp = dbValues[_MaxwellIsotropic3D::didVp];
  const double viscosity = dbValues[_MaxwellIsotropic3D::didViscosity];
 
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
  const double maxwelltime = viscosity / mu;

  paramVals[_MaxwellIsotropic3D::pidDensity] = density;
  paramVals[_MaxwellIsotropic3D::pidMu] = mu;
  paramVals[_MaxwellIsotropic3D::pidLambda] = lambda;
  paramVals[_MaxwellIsotropic3D::pidMaxwellTime] = maxwelltime;

  PetscLogFlopsNoCheck(7);
} // _dbToParameters

// ----------------------------------------------------------------------
// Get number of entries in stress tensor.
int
pylith::materials::MaxwellIsotropic3D::_tensorSize(void) const
{ // _tensorSize
  return _MaxwellIsotropic3D::tensorSize;
} // _tensorSize

// ----------------------------------------------------------------------
// Get number of entries in elasticity matrix for material.
int
pylith::materials::MaxwellIsotropic3D::_numElasticConsts(void) const
{ // numElasticConsts
  return _MaxwellIsotropic3D::numElasticConsts;
} // numElasticConsts

// ----------------------------------------------------------------------
// Compute density at location from parameters.
void
pylith::materials::MaxwellIsotropic3D::_calcDensity(double* const density,
						    const double* parameters,
						    const int numParams)
{ // _calcDensity
  assert(0 != density);
  assert(0 != parameters);
  assert(_numParamsQuadPt == numParams);

  density[0] = parameters[_MaxwellIsotropic3D::pidDensity];
} // _calcDensity

// ----------------------------------------------------------------------
// Compute stress tensor at location from parameters as an elastic
// material.
void
pylith::materials::MaxwellIsotropic3D::_calcStressElastic(
						    double* const stress,
						    const int stressSize,
						    const double* parameters,
						    const int numParams,
						    const double* totalStrain,
						    const int strainSize)
{ // _calcStressElastic
  assert(0 != stress);
  assert(_MaxwellIsotropic3D::tensorSize == stressSize);
  assert(0 != parameters);
  assert(_numParamsQuadPt == numParams);
  assert(0 != totalStrain);
  assert(_MaxwellIsotropic3D::tensorSize == strainSize);

  const double density = parameters[_MaxwellIsotropic3D::pidDensity];
  const double mu = parameters[_MaxwellIsotropic3D::pidMu];
  const double lambda = parameters[_MaxwellIsotropic3D::pidLambda];
  const double maxwelltime = parameters[_MaxwellIsotropic3D::pidMaxwellTime];

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
  // for (int iComp=0; iComp < _MaxwellIsotropic3D::tensorSize; ++iComp)
    // std::cout << "  " << totalStrain[iComp];
  // std::cout << std::endl << " stress: " << std::endl;
  // for (int iComp=0; iComp < _MaxwellIsotropic3D::tensorSize; ++iComp)
    // std::cout << "  " << (*stress)[iComp];
  // std::cout << std::endl;

  PetscLogFlopsNoCheck(13);
} // _calcStressElastic

// ----------------------------------------------------------------------
// Compute stress tensor at location from parameters as a viscoelastic
// material.
void
pylith::materials::MaxwellIsotropic3D::_calcStressViscoelastic(
						    double* const stress,
						    const int stressSize,
						    const double* parameters,
						    const int numParams,
						    const double* totalStrain,
						    const int strainSize)
{ // _calcStressElastic
  assert(0 != stress);
  assert(_MaxwellIsotropic3D::tensorSize == stressSize);
  assert(0 != parameters);
  assert(_numParamsQuadPt == numParams);
  assert(0 != totalStrain);
  assert(_MaxwellIsotropic3D::tensorSize == strainSize);

  const double density = parameters[_MaxwellIsotropic3D::pidDensity];
  const double mu = parameters[_MaxwellIsotropic3D::pidMu];
  const double lambda = parameters[_MaxwellIsotropic3D::pidLambda];
  const double maxwelltime = parameters[_MaxwellIsotropic3D::pidMaxwellTime];

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
  // See what's going on in state variables.
  // std::cout << " pidStrainT, pidVisStrain : " << std::endl;
  // for (int iComp=0; iComp < _MaxwellIsotropic3D::tensorSize; ++iComp)
    // std::cout << "  " << parameters[_MaxwellIsotropic3D::pidStrainT][iComp]
	    // << "   " << parameters[_MaxwellIsotropic3D::pidVisStrain][iComp]
	    // << std::endl;

  const double meanStrainT = (parameters[_MaxwellIsotropic3D::pidStrainT+0] +
			      parameters[_MaxwellIsotropic3D::pidStrainT+1] +
			      parameters[_MaxwellIsotropic3D::pidStrainT+2])/3.0;
  
  PetscLogFlopsNoCheck(11);

  // Time integration.
  double dq = ViscoelasticMaxwell::computeVisStrain(_dt, maxwelltime);
  const double expFac = exp(-_dt/maxwelltime);
  const double elasFac = 2.0*mu;
  double devStrainTpdt = 0.0;
  double devStrainT = 0.0;
  double devStressTpdt = 0.0;
  double visStrain = 0.0;
  PetscLogFlopsNoCheck(4);
  // std::cout << " _calcStressViscoelastic: " << std::endl;
  // std::cout << " stress  totalStrain  visStrain: " << std::endl;
  for (int iComp=0; iComp < _MaxwellIsotropic3D::tensorSize; ++iComp) {
    devStrainTpdt = totalStrain[iComp] - diag[iComp]*meanStrainTpdt;
    devStrainT = parameters[_MaxwellIsotropic3D::pidStrainT+iComp] -
      diag[iComp]*meanStrainT;
    visStrain = expFac*parameters[_MaxwellIsotropic3D::pidVisStrain+iComp] +
      dq*(devStrainTpdt - devStrainT);
    devStressTpdt = elasFac*visStrain;
    // Later I will want to put in initial stresses.
    stress[iComp] =diag[iComp]*meanStressTpdt+devStressTpdt;

    // Temporary to get stresses and strains.
    // std::cout << "  " << stress[iComp] << "  " << totalStrain[iComp] << "  " << visStrain << std:: endl;
  } // for

  PetscLogFlopsNoCheck(11 * _MaxwellIsotropic3D::tensorSize);
} // _calcStress

// ----------------------------------------------------------------------
// Compute derivative of elasticity matrix at location from parameters.
void
pylith::materials::MaxwellIsotropic3D::_calcElasticConstsElastic(
						  double* const elasticConsts,
						  const int numElasticConsts,
						  const double* parameters,
						  const int numParams,
						  const double* totalStrain,
						  const int strainSize)
{ // _calcElasticConstsElastic
  assert(0 != elasticConsts);
  assert(_MaxwellIsotropic3D::numElasticConsts == numElasticConsts);
  assert(0 != parameters);
  assert(_numParamsQuadPt == numParams);
  assert(0 != totalStrain);
  assert(_MaxwellIsotropic3D::tensorSize == strainSize);
 
  const double density = parameters[_MaxwellIsotropic3D::pidDensity];
  const double mu = parameters[_MaxwellIsotropic3D::pidMu];
  const double lambda = parameters[_MaxwellIsotropic3D::pidLambda];
  const double maxwelltime = parameters[_MaxwellIsotropic3D::pidMaxwellTime];

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
pylith::materials::MaxwellIsotropic3D::_calcElasticConstsViscoelastic(
						  double* const elasticConsts,
						  const int numElasticConsts,
						  const double* parameters,
						  const int numParams,
						  const double* totalStrain,
						  const int strainSize)
{ // _calcElasticConstsViscoelastic
  assert(0 != elasticConsts);
  assert(_MaxwellIsotropic3D::numElasticConsts == numElasticConsts);
  assert(0 != parameters);
  assert(_numParamsQuadPt == numParams);
  assert(0 != totalStrain);
  assert(_MaxwellIsotropic3D::tensorSize == strainSize);
 
  const double density = parameters[_MaxwellIsotropic3D::pidDensity];
  const double mu = parameters[_MaxwellIsotropic3D::pidMu];
  const double lambda = parameters[_MaxwellIsotropic3D::pidLambda];
  const double maxwelltime = parameters[_MaxwellIsotropic3D::pidMaxwellTime];

  const double mu2 = 2.0 * mu;
  const double bulkmodulus = lambda + mu2/3.0;

  double dq = ViscoelasticMaxwell::computeVisStrain(_dt, maxwelltime);

  const double visFac = mu*dq/3.0;
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

  PetscLogFlopsNoCheck(7);
} // _calcElasticConstsViscoelastic

// ----------------------------------------------------------------------
// Update state variables.
void
pylith::materials::MaxwellIsotropic3D::_updateStateElastic(
						 double* const parameters,
						 const int numParams,
						 const double* totalStrain,
						 const int strainSize)
{ // _updateStateElastic
  assert(0 != parameters);
  assert(_numParamsQuadPt == numParams);
  assert(0 != totalStrain);
  assert(_MaxwellIsotropic3D::tensorSize == strainSize);

  const double maxwelltime = parameters[_MaxwellIsotropic3D::pidMaxwellTime];

  const double e11 = totalStrain[0];
  const double e22 = totalStrain[1];
  const double e33 = totalStrain[2];

  const double traceStrainTpdt = e11 + e22 + e33;
  const double meanStrainTpdt = traceStrainTpdt/3.0;

  const double diag[] = { 1.0, 1.0, 1.0, 0.0, 0.0, 0.0 };

  // Temporary to get stresses.
  // double_array stress(6);
  // _calcStressElastic(&stress, (*parameters), totalStrain);

  for (int iComp=0; iComp < _MaxwellIsotropic3D::tensorSize; ++iComp) {
    parameters[_MaxwellIsotropic3D::pidStrainT+iComp] = totalStrain[iComp];
    parameters[_MaxwellIsotropic3D::pidVisStrain+iComp] =
      totalStrain[iComp] - diag[iComp]*meanStrainTpdt;
  } // for
  PetscLogFlopsNoCheck(5 * _MaxwellIsotropic3D::tensorSize);
//   std::cout << std::endl;
//   std::cout << " updateStateElastic: "<< std::endl;
//   std::cout << " StrainT  VisStrain  Stress: " << std::endl;
//   for (int iComp=0; iComp < _MaxwellIsotropic3D::tensorSize; ++iComp)
//     std::cout << "  " << parameters[_MaxwellIsotropic3D::pidStrainT+iComp]
// 	    << "   " << parameters[_MaxwellIsotropic3D::pidVisStrain+iComp]
// 	    << "   " << stress[iComp]
// 	    << std::endl;

  _needNewJacobian = true;
} // _updateStateElastic

// ----------------------------------------------------------------------
// Update state variables.
void
pylith::materials::MaxwellIsotropic3D::_updateStateViscoelastic(
						 double* const parameters,
						 const int numParams,
						 const double* totalStrain,
						 const int strainSize)
{ // _updateStateViscoelastic
  assert(0 != parameters);
  assert(_numParamsQuadPt == numParams);
  assert(0 != totalStrain);
  assert(_MaxwellIsotropic3D::tensorSize == strainSize);

  const double maxwelltime = parameters[_MaxwellIsotropic3D::pidMaxwellTime];

  const double e11 = totalStrain[0];
  const double e22 = totalStrain[1];
  const double e33 = totalStrain[2];

  const double meanStrainTpdt = (e11 + e22 + e33) / 3.0;

  const double diag[] = { 1.0, 1.0, 1.0, 0.0, 0.0, 0.0 };

  // Temporary to get stresses.
  double stress[6];
  const int stressSize = 6;
  _calcStressViscoelastic(stress, stressSize, 
			  parameters, numParams,
			  totalStrain, strainSize);

  const double meanStrainT = 
    (parameters[_MaxwellIsotropic3D::pidStrainT+0] +
     parameters[_MaxwellIsotropic3D::pidStrainT+1] +
     parameters[_MaxwellIsotropic3D::pidStrainT+2]) / 3.0;
  
  PetscLogFlopsNoCheck(6);

  double dq = ViscoelasticMaxwell::computeVisStrain(_dt, maxwelltime);

  const double expFac = exp(-_dt/maxwelltime);
  double devStrainTpdt = 0.0;
  double devStrainT = 0.0;
  double visStrain = 0.0;
  PetscLogFlopsNoCheck(3);
  for (int iComp=0; iComp < _MaxwellIsotropic3D::tensorSize; ++iComp) {
    devStrainTpdt = totalStrain[iComp] - diag[iComp]*meanStrainTpdt;
    devStrainT = parameters[_MaxwellIsotropic3D::pidStrainT+iComp] -
      diag[iComp] * meanStrainT;
    visStrain = expFac * 
      parameters[_MaxwellIsotropic3D::pidVisStrain+iComp] +
      dq * (devStrainTpdt - devStrainT);
    parameters[_MaxwellIsotropic3D::pidVisStrain+iComp] = visStrain;
    parameters[_MaxwellIsotropic3D::pidStrainT+iComp] = totalStrain[iComp];
  } // for
  PetscLogFlopsNoCheck(8 * _MaxwellIsotropic3D::tensorSize);

  _needNewJacobian = false;

//   std::cout << std::endl;
//   std::cout << " updateStateViscoelastic: "<< std::endl;
//   std::cout << " StrainT  VisStrain  Stress: " << std::endl;
//   for (int iComp=0; iComp < _MaxwellIsotropic3D::tensorSize; ++iComp)
//     std::cout << "  " << parameters[_MaxwellIsotropic3D::pidStrainT+iComp]
// 	    << "   " << parameters[_MaxwellIsotropic3D::pidVisStrain+iComp]
// 	    << "   " << stress[iComp]
// 	    << std::endl;
} // _updateStateViscoelastic


// End of file 
