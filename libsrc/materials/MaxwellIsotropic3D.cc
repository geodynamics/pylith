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

#include "pylith/utils/array.hh" // USES double_array

#include <assert.h> // USES assert()
#include <sstream> // USES std::ostringstream
#include <stdexcept> // USES std::runtime_error

// ----------------------------------------------------------------------
namespace pylith {
  namespace materials {
    namespace _MaxwellIsotropic3D{;

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
      const char* namesParameters[] =
        {"density", "mu", "lambda" , "maxwellTime", "strainT", "visStrain"};

      /// Indices (order) of parameters
      const int pidDensity = 0;
      const int pidMu = 1;
      const int pidLambda = 2;
      const int pidMaxwellTime = 3;
      const int pidStrainT = 4;
      const int pidVisStrain = 5;
    } // _MaxwellIsotropic3D
  } // materials
} // pylith

// ----------------------------------------------------------------------
// Default constructor.
pylith::materials::MaxwellIsotropic3D::MaxwellIsotropic3D(void)
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
// Get names of parameters for physical properties.
const char**
pylith::materials::MaxwellIsotropic3D::_parameterNames(void) const
{ // _parameterNames
  return _MaxwellIsotropic3D::namesParameters;
} // _parameterNames

// ----------------------------------------------------------------------
// Get number of values for each parameter
void
pylith::materials::MaxwellIsotropic3D::_numParamValues(int_array* numValues) const
{ // _numParamValues
  assert(0 != numValues);

  const int numParams = _MaxwellIsotropic3D::numParameters;
  numValues->resize(numParams);
  for (int i=0; i< numParams; ++i)
    (*numValues)[i] = _MaxwellIsotropic3D::numParamValues[i];
} // _numParamValues

// ----------------------------------------------------------------------
// Compute parameters from values in spatial database.
void
pylith::materials::MaxwellIsotropic3D::_dbToParameters(std::vector<double_array>* paramVals,
					  const double_array& dbValues) const
{ // _dbToParameters
  assert(0 != paramVals);
  const int numParams = paramVals->size();
  assert(_MaxwellIsotropic3D::numParameters == numParams);
  const int numDBValues = dbValues.size();
  assert(_MaxwellIsotropic3D::numDBValues == numDBValues);
  for (int i=0; i < numParams; ++i)
    assert(_MaxwellIsotropic3D::numParamValues[i] == (*paramVals)[i].size());

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

  (*paramVals)[_MaxwellIsotropic3D::pidDensity][0] = density;
  (*paramVals)[_MaxwellIsotropic3D::pidMu][0] = mu;
  (*paramVals)[_MaxwellIsotropic3D::pidLambda][0] = lambda;
  (*paramVals)[_MaxwellIsotropic3D::pidMaxwellTime][0] = maxwelltime;
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
pylith::materials::MaxwellIsotropic3D::_calcDensity(
				double_array* const density,
				const std::vector<double_array>& parameters)
{ // _calcDensity
  assert(0 != density);
  assert(1 == density->size());
  assert(_MaxwellIsotropic3D::numParameters == parameters.size());
  assert(1 == parameters[_MaxwellIsotropic3D::pidDensity].size());

  (*density)[0] = parameters[_MaxwellIsotropic3D::pidDensity][0];
} // _calcDensity

// ----------------------------------------------------------------------
// Compute stress tensor at location from parameters.
void
pylith::materials::MaxwellIsotropic3D::_calcStress(double_array* const stress,
						   const std::vector<double_array>& parameters,
						   const double_array& totalStrain)
{ // _calcStress
  assert(0 != stress);
  assert(_MaxwellIsotropic3D::tensorSize == stress->size());
  assert(_MaxwellIsotropic3D::numParameters == parameters.size());
  assert(_MaxwellIsotropic3D::tensorSize == totalStrain.size());
  assert(1 == parameters[_MaxwellIsotropic3D::pidDensity].size());
  assert(1 == parameters[_MaxwellIsotropic3D::pidMu].size());
  assert(1 == parameters[_MaxwellIsotropic3D::pidLambda].size());
  assert(6 == parameters[_MaxwellIsotropic3D::pidStrainT].size());
  assert(6 == parameters[_MaxwellIsotropic3D::pidVisStrain].size());

  const double density = parameters[_MaxwellIsotropic3D::pidDensity][0];
  const double mu = parameters[_MaxwellIsotropic3D::pidMu][0];
  const double lambda = parameters[_MaxwellIsotropic3D::pidLambda][0];
  const double maxwelltime = parameters[_MaxwellIsotropic3D::pidMaxwellTime][0];

  const double lambda2mu = lambda + 2.0 * mu;
  const double bulkmodulus = lambda + 2.0 * mu/3.0;

  const double e11 = totalStrain[0];
  const double e22 = totalStrain[1];
  const double e33 = totalStrain[2];
  const double e12 = totalStrain[3];
  const double e23 = totalStrain[4];
  const double e13 = totalStrain[5];
  
  const double traceStrainTpdt = e11 + e22 + e33;
  const double meanStrainTpdt = traceStrainTpdt/3.0;
  const double s123 = lambda * traceStrainTpdt;

  if (useElasticBehavior()) {
    (*stress)[0] = s123 + 2.0*mu*e11;
    (*stress)[1] = s123 + 2.0*mu*e22;
    (*stress)[2] = s123 + 2.0*mu*e33;
    (*stress)[3] = 2.0 * mu * e12;
    (*stress)[4] = 2.0 * mu * e23;
    (*stress)[5] = 2.0 * mu * e13;
  } else {
    const double diag[] = { 1.0, 1.0, 1.0, 0.0, 0.0, 0.0 };
    const double meanStressTpdt = bulkmodulus * traceStrainTpdt;
    const double meanStrainT = parameters[_MaxwellIsotropic3D::pidStrainT][0] +
      parameters[_MaxwellIsotropic3D::pidStrainT][1] +
      parameters[_MaxwellIsotropic3D::pidStrainT][2];

    // The code below should probably be in a separate function since it
    // is used more than once.  I should also probably cover the possibility
    // that Maxwell time is zero (although this should never happen).
    const double timeFrac = 1.0e-5;
    const int numTerms = 5;
    double dq = 0.0;
    if(maxwelltime < timeFrac*_dt) {
      double fSign = 1.0;
      double factorial = 1.0;
      double fraction = 1.0;
      dq = 1.0;
      for (int iTerm=2; iTerm <= numTerms; ++iTerm) {
	factorial *= iTerm;
	fSign *= -1.0;
	fraction *= _dt/maxwelltime;
	dq += fSign*fraction/factorial;
      } // for
    } else
      dq = maxwelltime*(1.0-exp(-_dt/maxwelltime))/_dt;
    const double expFac = exp(-_dt/maxwelltime);
    const double elasFac = 2.0*mu;
    double devStrainTpdt = 0.0;
    double devStrainT = 0.0;
    double devStressTpdt = 0.0;
    double visStrain = 0.0;
    for (int iComp=0; iComp < _MaxwellIsotropic3D::tensorSize; ++iComp) {
      devStrainTpdt = totalStrain[iComp] - diag[iComp]*meanStrainTpdt;
      devStrainT = parameters[_MaxwellIsotropic3D::pidStrainT][iComp] -
	diag[iComp]*meanStrainT;
      visStrain = expFac*parameters[_MaxwellIsotropic3D::pidVisStrain][iComp] +
	dq*(devStrainTpdt - devStrainT);
      devStressTpdt = elasFac*visStrain;
      // Later I will want to put in initial stresses.
      (*stress)[iComp] =diag[iComp]*meanStressTpdt+devStressTpdt;
    } // for
  } //else
} // _calcStress

// ----------------------------------------------------------------------
// Compute derivative of elasticity matrix at location from parameters.
void
pylith::materials::MaxwellIsotropic3D::_calcElasticConsts(
				       double_array* const elasticConsts,
				       const std::vector<double_array>& parameters,
				       const double_array& totalStrain)
{ // _calcElasticConsts
  assert(0 != elasticConsts);
  assert(_MaxwellIsotropic3D::numElasticConsts == elasticConsts->size());
  assert(_MaxwellIsotropic3D::numParameters == parameters.size());
  assert(_MaxwellIsotropic3D::tensorSize == totalStrain.size());
 
  const double density = parameters[_MaxwellIsotropic3D::pidDensity][0];
  const double mu = parameters[_MaxwellIsotropic3D::pidMu][0];
  const double lambda = parameters[_MaxwellIsotropic3D::pidLambda][0];
  const double maxwelltime = parameters[_MaxwellIsotropic3D::pidMaxwellTime][0];

  const double lambda2mu = lambda + 2.0 * mu;
  const double bulkmodulus = lambda + 2.0 * mu/3.0;

  if (useElasticBehavior()) {
    (*elasticConsts)[ 0] = lambda2mu; // C1111
    (*elasticConsts)[ 1] = lambda; // C1122
    (*elasticConsts)[ 2] = lambda; // C1133
    (*elasticConsts)[ 3] = 0; // C1112
    (*elasticConsts)[ 4] = 0; // C1123
    (*elasticConsts)[ 5] = 0; // C1113
    (*elasticConsts)[ 6] = lambda2mu; // C2222
    (*elasticConsts)[ 7] = lambda; // C2233
    (*elasticConsts)[ 8] = 0; // C2212
    (*elasticConsts)[ 9] = 0; // C2223
    (*elasticConsts)[10] = 0; // C2213
    (*elasticConsts)[11] = lambda2mu; // C3333
    (*elasticConsts)[12] = 0; // C3312
    (*elasticConsts)[13] = 0; // C3323
    (*elasticConsts)[14] = 0; // C3313
    (*elasticConsts)[15] = 2.0 * mu; // C1212
    (*elasticConsts)[16] = 0; // C1223
    (*elasticConsts)[17] = 0; // C1213
    (*elasticConsts)[18] = 2.0 * mu; // C2323
    (*elasticConsts)[19] = 0; // C2313
    (*elasticConsts)[20] = 2.0 * mu; // C1313
  } else {
    const double timeFrac = 1.0e-5;
    const int numTerms = 5;
    double dq = 0.0;
    if(maxwelltime < timeFrac*_dt) {
      double fSign = 1.0;
      double factorial = 1.0;
      double fraction = 1.0;
      dq = 1.0;
      for (int iTerm=2; iTerm <= numTerms; ++iTerm) {
	factorial *= iTerm;
	fSign *= -1.0;
	fraction *= _dt/maxwelltime;
	dq += fSign*fraction/factorial;
      } // for
    } else
      dq = maxwelltime*(1.0-exp(-_dt/maxwelltime))/_dt;
    const double visFac = mu*dq/3.0;
    (*elasticConsts)[ 0] = bulkmodulus + 4.0*dq; // C1111
    (*elasticConsts)[ 1] = bulkmodulus - 2.0*dq; // C1122
    (*elasticConsts)[ 2] = (*elasticConsts)[1]; // C1133
    (*elasticConsts)[ 3] = 0; // C1112
    (*elasticConsts)[ 4] = 0; // C1123
    (*elasticConsts)[ 5] = 0; // C1113
    (*elasticConsts)[ 6] = (*elasticConsts)[0]; // C2222
    (*elasticConsts)[ 7] = (*elasticConsts)[1]; // C2233
    (*elasticConsts)[ 8] = 0; // C2212
    (*elasticConsts)[ 9] = 0; // C2223
    (*elasticConsts)[10] = 0; // C2213
    (*elasticConsts)[11] = (*elasticConsts)[0]; // C3333
    (*elasticConsts)[12] = 0; // C3312
    (*elasticConsts)[13] = 0; // C3323
    (*elasticConsts)[14] = 0; // C3313
    (*elasticConsts)[15] = 3.0 * visFac; // C1212
    (*elasticConsts)[16] = 0; // C1223
    (*elasticConsts)[17] = 0; // C1213
    (*elasticConsts)[18] = (*elasticConsts)[15]; // C2323
    (*elasticConsts)[19] = 0; // C2313
    (*elasticConsts)[20] = (*elasticConsts)[15]; // C1313
  } // else
   
} // _calcElasticConsts

// ----------------------------------------------------------------------
// Update state variables.
void
pylith::materials::MaxwellIsotropic3D::_updateState(
				std::vector<double_array>& parameters,
				const double_array& totalStrain)
{ // _updateState
  assert(_MaxwellIsotropic3D::numParameters == parameters.size());
  assert(_MaxwellIsotropic3D::tensorSize == totalStrain.size());
  assert(1 == parameters[_MaxwellIsotropic3D::pidDensity].size());
  assert(1 == parameters[_MaxwellIsotropic3D::pidMu].size());
  assert(1 == parameters[_MaxwellIsotropic3D::pidLambda].size());
  assert(6 == parameters[_MaxwellIsotropic3D::pidStrainT].size());
  assert(6 == parameters[_MaxwellIsotropic3D::pidVisStrain].size());

  const double density = parameters[_MaxwellIsotropic3D::pidDensity][0];
  const double mu = parameters[_MaxwellIsotropic3D::pidMu][0];
  const double lambda = parameters[_MaxwellIsotropic3D::pidLambda][0];
  const double maxwelltime = parameters[_MaxwellIsotropic3D::pidMaxwellTime][0];

  const double lambda2mu = lambda + 2.0 * mu;
  const double bulkmodulus = lambda + 2.0 * mu/3.0;

  const double e11 = totalStrain[0];
  const double e22 = totalStrain[1];
  const double e33 = totalStrain[2];

  const double traceStrainTpdt = e11 + e22 + e33;
  const double meanStrainTpdt = traceStrainTpdt/3.0;
  const double s123 = lambda * traceStrainTpdt;

  const double diag[] = { 1.0, 1.0, 1.0, 0.0, 0.0, 0.0 };

  if (useElasticBehavior()) {
    for (int iComp=0; iComp < _MaxwellIsotropic3D::tensorSize; ++iComp) {
      parameters[_MaxwellIsotropic3D::pidStrainT][iComp] = totalStrain[iComp];
      parameters[_MaxwellIsotropic3D::pidVisStrain][iComp] =
	totalStrain[iComp] - diag[iComp]*meanStrainTpdt;
    } // for
    useElasticBehavior(false);
  } else {
    const double meanStrainT = parameters[_MaxwellIsotropic3D::pidStrainT][0] +
      parameters[_MaxwellIsotropic3D::pidStrainT][1] +
      parameters[_MaxwellIsotropic3D::pidStrainT][2];

    // The code below should probably be in a separate function since it
    // is used more than once.  I should also probably cover the possibility
    // that Maxwell time is zero (although this should never happen).
    const double timeFrac = 1.0e-5;
    const int numTerms = 5;
    double dq = 0.0;
    if(maxwelltime < timeFrac*_dt) {
      double fSign = 1.0;
      double factorial = 1.0;
      double fraction = 1.0;
      dq = 1.0;
      for (int iTerm=2; iTerm <= numTerms; ++iTerm) {
	factorial *= iTerm;
	fSign *= -1.0;
	fraction *= _dt/maxwelltime;
	dq += fSign*fraction/factorial;
      } // for
    } else
      dq = maxwelltime*(1.0-exp(-_dt/maxwelltime))/_dt;
    const double expFac = exp(-_dt/maxwelltime);
    double devStrainTpdt = 0.0;
    double devStrainT = 0.0;
    double visStrain = 0.0;
    for (int iComp=0; iComp < _MaxwellIsotropic3D::tensorSize; ++iComp) {
      devStrainTpdt = totalStrain[iComp] - diag[iComp]*meanStrainTpdt;
      devStrainT = parameters[_MaxwellIsotropic3D::pidStrainT][iComp] -
	diag[iComp]*meanStrainT;
      visStrain = expFac*parameters[_MaxwellIsotropic3D::pidVisStrain][iComp] +
	dq*(devStrainTpdt - devStrainT);
      parameters[_MaxwellIsotropic3D::pidVisStrain][iComp] = visStrain;
      parameters[_MaxwellIsotropic3D::pidStrainT][iComp] = totalStrain[iComp];
    } // for
  } //else
} // _calcStress


// End of file 
