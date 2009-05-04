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
#include "pylith/utils/constdefs.h" // USES MAXDOUBLE

#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include "petsc.h" // USES PetscLogFlops

#include <cassert> // USES assert()
#include <cstring> // USES memcpy()
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

      /// Number of physical properties.
      const int numProperties = 7;
      
      /// Physical properties.
      const Material::PropMetaData properties[] = {
	{ "density", 1, OTHER_FIELD },
	{ "mu", 1, OTHER_FIELD },
	{ "lambda", 1, OTHER_FIELD },
	{ "shear_ratio", numMaxwellModels, OTHER_FIELD },
	{ "maxwell_time", numMaxwellModels, OTHER_FIELD },
	{ "total_strain", tensorSize, OTHER_FIELD },
	{ "viscous_strain", numMaxwellModels*tensorSize, OTHER_FIELD },
      };
      /// Indices (order) of properties.
      const int pidDensity = 0;
      const int pidMuTot = pidDensity + 1;
      const int pidLambdaTot = pidMuTot + 1;
      const int pidShearRatio = pidLambdaTot + 1;
      const int pidMaxwellTime = pidShearRatio + numMaxwellModels;
      const int pidStrainT = pidMaxwellTime + numMaxwellModels;
      const int pidVisStrain = pidStrainT + tensorSize;      

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

      /// Indices (order) of database values.
      const int didDensity = 0;
      const int didVs = 1;
      const int didVp = 2;
      const int didShearRatio1 = 3;
      const int didShearRatio2 = 4;
      const int didShearRatio3 = 5;
      const int didViscosity1 = 6;
      const int didViscosity2 = 7;
      const int didViscosity3 = 8;

      /// Initial state values expected in spatial database
      const int numInitialStateDBValues = tensorSize;
      const char* namesInitialStateDBValues[] = { "stress_xx", "stress_yy",
                                                  "stress_zz", "stress_xy",
                                                  "stress_yz", "stress_xz" };

    } // _GenMaxwellIsotropic3D
  } // materials
} // pylith

// ----------------------------------------------------------------------
// Default constructor.
pylith::materials::GenMaxwellIsotropic3D::GenMaxwellIsotropic3D(void) :
  ElasticMaterial(_GenMaxwellIsotropic3D::tensorSize,
		  _GenMaxwellIsotropic3D::numElasticConsts,
		  _GenMaxwellIsotropic3D::namesDBValues,
		  _GenMaxwellIsotropic3D::namesInitialStateDBValues,
		  _GenMaxwellIsotropic3D::numDBValues,
		  _GenMaxwellIsotropic3D::properties,
		  _GenMaxwellIsotropic3D::numProperties),
  _calcElasticConstsFn(&pylith::materials::GenMaxwellIsotropic3D::_calcElasticConstsElastic),
  _calcStressFn(&pylith::materials::GenMaxwellIsotropic3D::_calcStressElastic),
  _updatePropertiesFn(&pylith::materials::GenMaxwellIsotropic3D::_updatePropertiesElastic)
{ // constructor
  _dimension = 3;
  _visStrain.resize(_GenMaxwellIsotropic3D::numMaxwellModels *
		    _GenMaxwellIsotropic3D::tensorSize);
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::materials::GenMaxwellIsotropic3D::~GenMaxwellIsotropic3D(void)
{ // destructor
} // destructor

// ----------------------------------------------------------------------
// Compute parameters from values in spatial database.
void
pylith::materials::GenMaxwellIsotropic3D::_dbToProperties(
					    double* const propValues,
					    const double_array& dbValues) const
{ // _dbToProperties
  assert(0 != propValues);
  const int numDBValues = dbValues.size();
  assert(_GenMaxwellIsotropic3D::numDBValues == numDBValues);

  const double density = dbValues[_GenMaxwellIsotropic3D::didDensity];
  const double vs = dbValues[_GenMaxwellIsotropic3D::didVs];
  const double vp = dbValues[_GenMaxwellIsotropic3D::didVp];
  const int numMaxwellModels = _GenMaxwellIsotropic3D::numMaxwellModels;
 
  if (density <= 0.0 || vs <= 0.0 || vp <= 0.0) {
    std::ostringstream msg;
    msg << "Spatial database returned nonpositive value for physical "
	<< "properties.\n"
	<< "density: " << density << "\n"
	<< "vp: " << vp << "\n"
	<< "vs: " << vs << "\n";
    throw std::runtime_error(msg.str());
  } // if

  const double mu = density * vs*vs;
  const double lambda = density * vp*vp - 2.0*mu;

  if (lambda <= 0.0) {
    std::ostringstream msg;
    msg << "Attempted to set Lame's constant lambda to nonpositive value.\n"
	<< "density: " << density << "\n"
	<< "vp: " << vp << "\n"
	<< "vs: " << vs << "\n";
    throw std::runtime_error(msg.str());
  } // if

  assert(mu > 0);

  propValues[_GenMaxwellIsotropic3D::pidDensity] = density;
  propValues[_GenMaxwellIsotropic3D::pidMuTot] = mu;
  propValues[_GenMaxwellIsotropic3D::pidLambdaTot] = lambda;

  // Loop over number of Maxwell models.
  for (int model =0;
       model < numMaxwellModels; ++model) {
    double muRatio = dbValues[_GenMaxwellIsotropic3D::didShearRatio1 + model];
    double viscosity = dbValues[_GenMaxwellIsotropic3D::didViscosity1 + model];
    double muFac = muRatio*mu;
    double maxwellTime = 1.0e30;
    if (muFac > 0.0)
      maxwellTime = viscosity / muFac;
    if (muRatio < 0.0 || viscosity < 0.0 || muFac < 0.0 || maxwellTime < 0.0) {
      std::ostringstream msg;
      msg << "Found negative value(s) for physical properties.\n"
	  << "muRatio: " << muRatio << "\n"
	  << "viscosity: " << viscosity << "\n"
	  << "muFac: " << muFac << "\n"
	  << "maxwellTime: " << maxwellTime << "\n";
      throw std::runtime_error(msg.str());
    } // if
    propValues[_GenMaxwellIsotropic3D::pidShearRatio + model] = muRatio;
    propValues[_GenMaxwellIsotropic3D::pidMaxwellTime + model] = maxwellTime;
  } // for

  PetscLogFlops(6+2*numMaxwellModels);
} // _dbToProperties

// ----------------------------------------------------------------------
// Nondimensionalize properties.
void
pylith::materials::GenMaxwellIsotropic3D::_nondimProperties(double* const values,
							 const int nvalues) const
{ // _nondimProperties
  assert(0 != _normalizer);
  assert(0 != values);
  assert(nvalues == _totalPropsQuadPt);

  const double densityScale = _normalizer->densityScale();
  const double pressureScale = _normalizer->pressureScale();
  const double timeScale = _normalizer->timeScale();
  values[_GenMaxwellIsotropic3D::pidDensity] = 
    _normalizer->nondimensionalize(values[_GenMaxwellIsotropic3D::pidDensity],
				   densityScale);
  values[_GenMaxwellIsotropic3D::pidMuTot] = 
    _normalizer->nondimensionalize(values[_GenMaxwellIsotropic3D::pidMuTot],
				   pressureScale);
  values[_GenMaxwellIsotropic3D::pidLambdaTot] = 
    _normalizer->nondimensionalize(values[_GenMaxwellIsotropic3D::pidLambdaTot],
				   pressureScale);
  const int numMaxwellModels = _GenMaxwellIsotropic3D::numMaxwellModels;
  _normalizer->nondimensionalize(&values[_GenMaxwellIsotropic3D::pidMaxwellTime],
				 numMaxwellModels, timeScale);

  PetscLogFlops(3+1*numMaxwellModels);
} // _nondimProperties

// ----------------------------------------------------------------------
// Dimensionalize properties.
void
pylith::materials::GenMaxwellIsotropic3D::_dimProperties(double* const values,
						      const int nvalues) const
{ // _dimProperties
  assert(0 != _normalizer);
  assert(0 != values);
  assert(nvalues == _totalPropsQuadPt);

  const double densityScale = _normalizer->densityScale();
  const double pressureScale = _normalizer->pressureScale();
  const double timeScale = _normalizer->timeScale();
  values[_GenMaxwellIsotropic3D::pidDensity] = 
    _normalizer->dimensionalize(values[_GenMaxwellIsotropic3D::pidDensity],
				densityScale);
  values[_GenMaxwellIsotropic3D::pidMuTot] = 
    _normalizer->dimensionalize(values[_GenMaxwellIsotropic3D::pidMuTot],
				pressureScale);
  values[_GenMaxwellIsotropic3D::pidLambdaTot] = 
    _normalizer->dimensionalize(values[_GenMaxwellIsotropic3D::pidLambdaTot],
				pressureScale);
  const int numMaxwellModels = _GenMaxwellIsotropic3D::numMaxwellModels;
  _normalizer->dimensionalize(&values[_GenMaxwellIsotropic3D::pidMaxwellTime],
			      numMaxwellModels, timeScale);

  PetscLogFlops(3+1*numMaxwellModels);
} // _dimProperties

// ----------------------------------------------------------------------
// Nondimensionalize initial state.
void
pylith::materials::GenMaxwellIsotropic3D::_nondimInitState(double* const values,
							const int nvalues) const
{ // _nondimInitState
  assert(0 != _normalizer);
  assert(0 != values);
  assert(nvalues == _GenMaxwellIsotropic3D::numInitialStateDBValues);

  const double pressureScale = _normalizer->pressureScale();
  _normalizer->nondimensionalize(values, nvalues, pressureScale);

  PetscLogFlops(nvalues);
} // _nondimInitState

// ----------------------------------------------------------------------
// Dimensionalize initial state.
void
pylith::materials::GenMaxwellIsotropic3D::_dimInitState(double* const values,
						     const int nvalues) const
{ // _dimInitState
  assert(0 != _normalizer);
  assert(0 != values);
  assert(nvalues == _GenMaxwellIsotropic3D::numInitialStateDBValues);
  
  const double pressureScale = _normalizer->pressureScale();
  _normalizer->dimensionalize(values, nvalues, pressureScale);

  PetscLogFlops(nvalues);
} // _dimInitState

// ----------------------------------------------------------------------
// Compute density at location from properties.
void
pylith::materials::GenMaxwellIsotropic3D::_calcDensity(
					    double* const density,
					    const double* properties,
					    const int numProperties)
{ // _calcDensity
  assert(0 != density);
  assert(0 != properties);
  assert(_totalPropsQuadPt == numProperties);
  
  density[0] = properties[_GenMaxwellIsotropic3D::pidDensity];
} // _calcDensity

// ----------------------------------------------------------------------
// Compute viscous strain for current time step.
void
pylith::materials::GenMaxwellIsotropic3D::_computeStateVars(
					    const double* properties,
					    const int numProperties,
					    const double* totalStrain,
					    const int strainSize,
					    const double* initialState,
					    const int initialStateSize)
{ // _computeStateVars
  assert(0 != properties);
  assert(_totalPropsQuadPt == numProperties);
  assert(0 != totalStrain);
  assert(_GenMaxwellIsotropic3D::tensorSize == strainSize);
  assert(_GenMaxwellIsotropic3D::tensorSize == initialStateSize);

  const int tensorSize = _GenMaxwellIsotropic3D::tensorSize;
  const int numMaxwellModels = _GenMaxwellIsotropic3D::numMaxwellModels;

  const double muRatio[] = {
    properties[_GenMaxwellIsotropic3D::pidShearRatio],
    properties[_GenMaxwellIsotropic3D::pidShearRatio+1],
    properties[_GenMaxwellIsotropic3D::pidShearRatio+2]};
  const double maxwellTime[] = {
    properties[_GenMaxwellIsotropic3D::pidMaxwellTime],
    properties[_GenMaxwellIsotropic3D::pidMaxwellTime+1],
    properties[_GenMaxwellIsotropic3D::pidMaxwellTime+2]};

  const double e11 = totalStrain[0];
  const double e22 = totalStrain[1];
  const double e33 = totalStrain[2];
  const double e12 = totalStrain[3];
  const double e23 = totalStrain[4];
  const double e13 = totalStrain[5];
  
  const double meanStrainTpdt = (e11 + e22 + e33)/3.0;

  const double diag[] = { 1.0, 1.0, 1.0, 0.0, 0.0, 0.0 };

  const double meanStrainT = 
    (properties[_GenMaxwellIsotropic3D::pidStrainT+0] +
     properties[_GenMaxwellIsotropic3D::pidStrainT+1] +
     properties[_GenMaxwellIsotropic3D::pidStrainT+2])/3.0;
  
  PetscLogFlops(6);

  // Compute Prony series terms
  double dq[] = { 0.0, 0.0, 0.0};
  for (int model=0; model < numMaxwellModels; ++model) {
    if (muRatio[model] != 0.0)
      dq[model] =
	ViscoelasticMaxwell::computeVisStrain(_dt, maxwellTime[model]);
  } // for

  // Compute new viscous strains
  double devStrainTpdt = 0.0;
  double devStrainT = 0.0;
  double deltaStrain = 0.0;
  
  for (int iComp=0; iComp < tensorSize; ++iComp) {
    devStrainTpdt = totalStrain[iComp] - diag[iComp]*meanStrainTpdt;
    devStrainT = properties[_GenMaxwellIsotropic3D::pidStrainT+iComp] -
      diag[iComp]*meanStrainT;
    deltaStrain = devStrainTpdt - devStrainT;
    PetscLogFlops(5);
    for (int model=0; model < numMaxwellModels; ++model) {
      if (muRatio[model] != 0.0) {
	int pindex = iComp + model * tensorSize;
	_visStrain[pindex] = exp(-_dt/maxwellTime[model])*
	  properties[_GenMaxwellIsotropic3D::pidVisStrain + pindex] +
	  dq[model] * deltaStrain;
	PetscLogFlops(6);
      } // if
    } // for
  } // for
} // _computeStateVars

// ----------------------------------------------------------------------
// Compute stress tensor at location from properties as an elastic
// material.
void
pylith::materials::GenMaxwellIsotropic3D::_calcStressElastic(
					    double* const stress,
					    const int stressSize,
					    const double* properties,
					    const int numProperties,
					    const double* totalStrain,
					    const int strainSize,
					    const double* initialState,
					    const int initialStateSize,
					    const bool computeStateVars)
{ // _calcStressElastic
  assert(0 != stress);
  assert(_GenMaxwellIsotropic3D::tensorSize == stressSize);
  assert(0 != properties);
  assert(_totalPropsQuadPt == numProperties);
  assert(0 != totalStrain);
  assert(_GenMaxwellIsotropic3D::tensorSize == strainSize);
  assert(_GenMaxwellIsotropic3D::tensorSize == initialStateSize);

  const double mu = properties[_GenMaxwellIsotropic3D::pidMuTot];
  const double lambda = properties[_GenMaxwellIsotropic3D::pidLambdaTot];
  const double mu2 = 2.0 * mu;

  const double e11 = totalStrain[0];
  const double e22 = totalStrain[1];
  const double e33 = totalStrain[2];
  const double e12 = totalStrain[3];
  const double e23 = totalStrain[4];
  const double e13 = totalStrain[5];
  
  const double traceStrainTpdt = e11 + e22 + e33;
  const double s123 = lambda * traceStrainTpdt;

  stress[0] = s123 + mu2*e11 + initialState[0];
  stress[1] = s123 + mu2*e22 + initialState[1];
  stress[2] = s123 + mu2*e33 + initialState[2];
  stress[3] = mu2 * e12 + initialState[3];
  stress[4] = mu2 * e23 + initialState[4];
  stress[5] = mu2 * e13 + initialState[5];
  // std::cout << " _calcStressElastic: " << std::endl;
  // std::cout << " totalStrain: " << std::endl;
  // for (int iComp=0; iComp < _GenMaxwellIsotropic3D::tensorSize; ++iComp)
    // std::cout << "  " << totalStrain[iComp];
  // std::cout << std::endl << " stress: " << std::endl;
  // for (int iComp=0; iComp < _GenMaxwellIsotropic3D::tensorSize; ++iComp)
    // std::cout << "  " << (*stress)[iComp];
  // std::cout << std::endl;

  PetscLogFlops(19);
} // _calcStressElastic


// ----------------------------------------------------------------------
// Compute stress tensor at location from properties as a viscoelastic
// material.
void
pylith::materials::GenMaxwellIsotropic3D::_calcStressViscoelastic(
					    double* const stress,
					    const int stressSize,
					    const double* properties,
					    const int numProperties,
					    const double* totalStrain,
					    const int strainSize,
					    const double* initialState,
					    const int initialStateSize,
					    const bool computeStateVars)
{ // _calcStressViscoelastic
  assert(0 != stress);
  assert(_GenMaxwellIsotropic3D::tensorSize == stressSize);
  assert(0 != properties);
  assert(_totalPropsQuadPt == numProperties);
  assert(0 != totalStrain);
  assert(_GenMaxwellIsotropic3D::tensorSize == strainSize);
  assert(_GenMaxwellIsotropic3D::tensorSize == initialStateSize);

  const int tensorSize = _GenMaxwellIsotropic3D::tensorSize;
  const int numMaxwellModels = _GenMaxwellIsotropic3D::numMaxwellModels;

  const double mu = properties[_GenMaxwellIsotropic3D::pidMuTot];
  const double lambda = properties[_GenMaxwellIsotropic3D::pidLambdaTot];
  const double muRatio[] = {
    properties[_GenMaxwellIsotropic3D::pidShearRatio],
    properties[_GenMaxwellIsotropic3D::pidShearRatio+1],
    properties[_GenMaxwellIsotropic3D::pidShearRatio+2]};

  const double mu2 = 2.0 * mu;
  const double bulkModulus = lambda + mu2/3.0;

  const double e11 = totalStrain[0];
  const double e22 = totalStrain[1];
  const double e33 = totalStrain[2];
  const double e12 = totalStrain[3];
  const double e23 = totalStrain[4];
  const double e13 = totalStrain[5];
  
  const double traceStrainTpdt = e11 + e22 + e33;
  const double meanStrainTpdt = traceStrainTpdt/3.0;
  const double meanStressTpdt = bulkModulus * traceStrainTpdt;

  const double diag[] = { 1.0, 1.0, 1.0, 0.0, 0.0, 0.0 };
  
  double visFrac = 0.0;
  for (int model = 0; model < numMaxwellModels; ++model) 
    visFrac += muRatio[model];

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

  PetscLogFlops(8 + numMaxwellModels);

  // Get viscous strains
  if (computeStateVars) {
    pylith::materials::GenMaxwellIsotropic3D::_computeStateVars(properties,
								numProperties,
								totalStrain,
								strainSize,
								initialState,
								initialStateSize);
  } else {
    memcpy(&_visStrain[0], &properties[_GenMaxwellIsotropic3D::pidVisStrain],
	   numMaxwellModels * tensorSize * sizeof(double));
  } // else


  // Compute new stresses
  double devStrainTpdt = 0.0;
  double devStressTpdt = 0.0;
  // std::cout << " _calcStressViscoelastic: " << std::endl;
  // std::cout << " stress  totalStrain  _visStrain: " << std::endl;
  for (int iComp=0; iComp < tensorSize; ++iComp) {
    devStrainTpdt = totalStrain[iComp] - diag[iComp]*meanStrainTpdt;
    devStressTpdt = elasFrac*devStrainTpdt;
    for (int model=0; model < numMaxwellModels; ++model) {
      if (muRatio[model] != 0.0) {
	int pindex = iComp + model * tensorSize;
	devStressTpdt += muRatio[model] * _visStrain[pindex];
      } // if
    } // for

    devStressTpdt = mu2*devStressTpdt;
    stress[iComp] =diag[iComp] * meanStressTpdt + devStressTpdt +
	    initialState[iComp];
    // std::cout << devStressTpdt << "  " << stress[iComp] << std::endl;

    // Temporary to get stresses and strains.
    // std::cout << "  " << stress[iComp] << "  " << totalStrain[iComp] << "  " << _visStrain << std:: endl;
  } // for

  PetscLogFlops((7 + 2*numMaxwellModels) * tensorSize);
} // _calcStressViscoelastic

// ----------------------------------------------------------------------
// Compute derivative of elasticity matrix at location from properties.
void
pylith::materials::GenMaxwellIsotropic3D::_calcElasticConstsElastic(
					    double* const elasticConsts,
					    const int numElasticConsts,
					    const double* properties,
					    const int numProperties,
					    const double* totalStrain,
					    const int strainSize,
					    const double* initialState,
					    const int initialStateSize)
{ // _calcElasticConstsElastic
  assert(0 != elasticConsts);
  assert(_GenMaxwellIsotropic3D::numElasticConsts == numElasticConsts);
  assert(0 != properties);
  assert(_totalPropsQuadPt == numProperties);
  assert(0 != totalStrain);
  assert(_GenMaxwellIsotropic3D::tensorSize == strainSize);
  assert(_GenMaxwellIsotropic3D::tensorSize == initialStateSize);
 
  const double mu = properties[_GenMaxwellIsotropic3D::pidMuTot];
  const double lambda = properties[_GenMaxwellIsotropic3D::pidLambdaTot];

  const double mu2 = 2.0 * mu;
  const double lambda2mu = lambda + mu2;

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

  PetscLogFlops(4);
} // _calcElasticConstsElastic

// ----------------------------------------------------------------------
// Compute derivative of elasticity matrix at location from properties
// as a viscoelastic material.
void
pylith::materials::GenMaxwellIsotropic3D::_calcElasticConstsViscoelastic(
					    double* const elasticConsts,
					    const int numElasticConsts,
					    const double* properties,
					    const int numProperties,
					    const double* totalStrain,
					    const int strainSize,
					    const double* initialState,
					    const int initialStateSize)
{ // _calcElasticConstsViscoelastic
  assert(0 != elasticConsts);
  assert(_GenMaxwellIsotropic3D::numElasticConsts == numElasticConsts);
  assert(0 != properties);
  assert(_totalPropsQuadPt == numProperties);
  assert(0 != totalStrain);
  assert(_GenMaxwellIsotropic3D::tensorSize == strainSize);
  assert(_GenMaxwellIsotropic3D::tensorSize == initialStateSize);
 
  const double mu = properties[_GenMaxwellIsotropic3D::pidMuTot];
  const double lambda = properties[_GenMaxwellIsotropic3D::pidLambdaTot];
  const int numMaxwellModels = _GenMaxwellIsotropic3D::numMaxwellModels;

  const double mu2 = 2.0 * mu;
  const double bulkModulus = lambda + mu2/3.0;

  // Compute viscous contribution.
  double visFac = 0.0;
  double visFrac = 0.0;
  double shearRatio = 0.0;
  for (int model = 0; model < numMaxwellModels; ++model) {
    shearRatio = properties[_GenMaxwellIsotropic3D::pidShearRatio + model];
    double maxwellTime = 1.0e30;
    visFrac += shearRatio;
    if (shearRatio != 0.0) {
      maxwellTime = properties[_GenMaxwellIsotropic3D::pidMaxwellTime + model];
      visFac +=
	shearRatio*ViscoelasticMaxwell::computeVisStrain(_dt, maxwellTime);
    } // if
  } // for
  double elasFrac = 1.0 - visFrac;
  double shearFac = mu*(elasFrac + visFac)/3.0;

  elasticConsts[ 0] = bulkModulus + 4.0*shearFac; // C1111
  elasticConsts[ 1] = bulkModulus - 2.0*shearFac; // C1122
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
  elasticConsts[15] = 6.0 * shearFac; // C1212
  elasticConsts[16] = 0; // C1223
  elasticConsts[17] = 0; // C1213
  elasticConsts[18] = elasticConsts[15]; // C2323
  elasticConsts[19] = 0; // C2313
  elasticConsts[20] = elasticConsts[15]; // C1313

#if 0 // DEBUGGING
  std::cout << "_calcElasticConstsViscoelastic" << std::endl;
  std::cout << elasticConsts[0] << "  " << elasticConsts[1] << "  " << elasticConsts[2] << std::endl;
  std::cout << elasticConsts[6] << "  " << elasticConsts[7] << std::endl;
  std::cout << elasticConsts[11] << std::endl;
  std::cout << elasticConsts[15] << std::endl;
  std::cout << elasticConsts[18] << std::endl;
  std::cout << elasticConsts[20] << std::endl;
#endif

  PetscLogFlops(8 + 2 * numMaxwellModels);
} // _calcElasticConstsViscoelastic

// ----------------------------------------------------------------------
// Get stable time step for implicit time integration.
double
pylith::materials::GenMaxwellIsotropic3D::_stableTimeStepImplicit(const double* properties,
				 const int numProperties) const
{ // _stableTimeStepImplicit
  assert(0 != properties);
  assert(_totalPropsQuadPt == numProperties);

  double dtStable = pylith::PYLITH_MAXDOUBLE;

  const int numMaxwellModels = _GenMaxwellIsotropic3D::numMaxwellModels;
  for (int i=0; i < numMaxwellModels; ++i) {
    const double maxwellTime = 
      properties[_GenMaxwellIsotropic3D::pidMaxwellTime+i];
    const double dt = 0.1*maxwellTime;
    if (dt < dtStable)
      dtStable = dt;
  } // for

  return dtStable;
} // _stableTimeStepImplicit

// ----------------------------------------------------------------------
// Update state variables.
void
pylith::materials::GenMaxwellIsotropic3D::_updatePropertiesElastic(
					    double* const properties,
					    const int numProperties,
					    const double* totalStrain,
					    const int strainSize,
					    const double* initialState,
					    const int initialStateSize)
{ // _updatePropertiesElastic
  assert(0 != properties);
  assert(_totalPropsQuadPt == numProperties);
  assert(0 != totalStrain);
  assert(_GenMaxwellIsotropic3D::tensorSize == strainSize);
  assert(_GenMaxwellIsotropic3D::tensorSize == initialStateSize);

  const int numMaxwellModels = _GenMaxwellIsotropic3D::numMaxwellModels;

  const double e11 = totalStrain[0];
  const double e22 = totalStrain[1];
  const double e33 = totalStrain[2];

  const double traceStrainTpdt = e11 + e22 + e33;
  const double meanStrainTpdt = traceStrainTpdt/3.0;

  const double diag[] = { 1.0, 1.0, 1.0, 0.0, 0.0, 0.0 };

  // Initialize all viscous strains to deviatoric elastic strains.
  double devStrain = 0.0;
  double shearRatio = 0.0;

  for (int iComp=0; iComp < _GenMaxwellIsotropic3D::tensorSize; ++iComp) {
    properties[_GenMaxwellIsotropic3D::pidStrainT+iComp] = totalStrain[iComp];
    devStrain = totalStrain[iComp] - diag[iComp] * meanStrainTpdt;
    for (int model = 0; model < numMaxwellModels; ++model) {
      shearRatio = properties[_GenMaxwellIsotropic3D::pidShearRatio + model];
      properties[_GenMaxwellIsotropic3D::pidVisStrain+
		 iComp+model*_GenMaxwellIsotropic3D::tensorSize] =
	devStrain;
    } // for
  } // for
  PetscLogFlops(3 + 2 * _GenMaxwellIsotropic3D::tensorSize);

  _needNewJacobian = true;
} // _updatePropertiesElastic

// ----------------------------------------------------------------------
// Update state variables.
void
pylith::materials::GenMaxwellIsotropic3D::_updatePropertiesViscoelastic(
					    double* const properties,
					    const int numProperties,
					    const double* totalStrain,
					    const int strainSize,
					    const double* initialState,
					    const int initialStateSize)
{ // _updatePropertiesViscoelastic
  assert(0 != properties);
  assert(_totalPropsQuadPt == numProperties);
  assert(0 != totalStrain);
  assert(_GenMaxwellIsotropic3D::tensorSize == strainSize);
  assert(_GenMaxwellIsotropic3D::tensorSize == initialStateSize);

  const int numMaxwellModels = _GenMaxwellIsotropic3D::numMaxwellModels;
  const int tensorSize = _GenMaxwellIsotropic3D::tensorSize;

  pylith::materials::GenMaxwellIsotropic3D::_computeStateVars(properties,
							      numProperties,
							      totalStrain,
							      strainSize,
							      initialState,
							      initialStateSize);

  memcpy(&properties[_GenMaxwellIsotropic3D::pidVisStrain],
	 &_visStrain[0], 
	 numMaxwellModels * tensorSize * sizeof(double));
  memcpy(&properties[_GenMaxwellIsotropic3D::pidStrainT],
	 &totalStrain[0], 
	 tensorSize * sizeof(double));

  _needNewJacobian = false;

} // _updatePropertiesViscoelastic


// End of file 
