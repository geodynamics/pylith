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

      /// Number of physical properties.
      const int numProperties = 6;

      /// Physical properties.
      const Material::PropMetaData properties[] = {
	{ "density", 1, SCALAR_FIELD },
	{ "lambda", 1, SCALAR_FIELD },
	{ "mu", 1, SCALAR_FIELD },
	{ "Maxwell-time", 1, SCALAR_FIELD },
	{ "total-strain", 6, TENSOR_FIELD },
	{ "viscous-strain", 6, TENSOR_FIELD },
      };
      /// Indices (order) of properties.
      const int pidDensity = 0;
      const int pidMu = pidDensity + 1;
      const int pidLambda = pidMu + 1;
      const int pidMaxwellTime = pidLambda + 1;
      const int pidStrainT = pidMaxwellTime + 1;
      const int pidVisStrain = pidStrainT + 6;

      /// Values expected in spatial database
      const int numDBValues = 4;
      const char* namesDBValues[] = {"density", "vs", "vp" , "viscosity"};

      /// Indices (order) of database values
      const int didDensity = 0;
      const int didVs = 1;
      const int didVp = 2;
      const int didViscosity = 3;

    } // _MaxwellIsotropic3D
  } // materials
} // pylith

// ----------------------------------------------------------------------
// Default constructor.
pylith::materials::MaxwellIsotropic3D::MaxwellIsotropic3D(void) :
  ElasticMaterial(_MaxwellIsotropic3D::tensorSize,
		  _MaxwellIsotropic3D::numElasticConsts,
		  _MaxwellIsotropic3D::namesDBValues,
		  _MaxwellIsotropic3D::numDBValues,
		  _MaxwellIsotropic3D::properties,
		  _MaxwellIsotropic3D::numProperties),
  _calcElasticConstsFn(&pylith::materials::MaxwellIsotropic3D::_calcElasticConstsElastic),
  _calcStressFn(&pylith::materials::MaxwellIsotropic3D::_calcStressElastic),
  _updatePropertiesFn(&pylith::materials::MaxwellIsotropic3D::_updatePropertiesElastic)
{ // constructor
  _dimension = 3;
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::materials::MaxwellIsotropic3D::~MaxwellIsotropic3D(void)
{ // destructor
} // destructor

// ----------------------------------------------------------------------
// Compute properties from values in spatial database.
void
pylith::materials::MaxwellIsotropic3D::_dbToProperties(
					    double* const propValues,
					    const double_array& dbValues) const
{ // _dbToProperties
  assert(0 != propValues);
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
  assert(maxwelltime > 0.0);

  propValues[_MaxwellIsotropic3D::pidDensity] = density;
  propValues[_MaxwellIsotropic3D::pidMu] = mu;
  propValues[_MaxwellIsotropic3D::pidLambda] = lambda;
  propValues[_MaxwellIsotropic3D::pidMaxwellTime] = maxwelltime;

  PetscLogFlopsNoCheck(7);
} // _dbToProperties

// ----------------------------------------------------------------------
// Compute density at location from properties.
void
pylith::materials::MaxwellIsotropic3D::_calcDensity(double* const density,
						    const double* properties,
						    const int numProperties)
{ // _calcDensity
  assert(0 != density);
  assert(0 != properties);
  assert(_totalPropsQuadPt == numProperties);

  density[0] = properties[_MaxwellIsotropic3D::pidDensity];
} // _calcDensity

// ----------------------------------------------------------------------
// Compute stress tensor at location from properties as an elastic
// material.
void
pylith::materials::MaxwellIsotropic3D::_calcStressElastic(
						  double* const stress,
						  const int stressSize,
						  const double* properties,
						  const int numProperties,
						  const double* totalStrain,
						  const int strainSize,
						  const bool computeStateVars)
{ // _calcStressElastic
  assert(0 != stress);
  assert(_MaxwellIsotropic3D::tensorSize == stressSize);
  assert(0 != properties);
  assert(_totalPropsQuadPt == numProperties);
  assert(0 != totalStrain);
  assert(_MaxwellIsotropic3D::tensorSize == strainSize);

  const double density = properties[_MaxwellIsotropic3D::pidDensity];
  const double mu = properties[_MaxwellIsotropic3D::pidMu];
  const double lambda = properties[_MaxwellIsotropic3D::pidLambda];
  const double maxwelltime = properties[_MaxwellIsotropic3D::pidMaxwellTime];

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
// Compute stress tensor at location from properties as a viscoelastic
// material.
void
pylith::materials::MaxwellIsotropic3D::_calcStressViscoelastic(
						  double* const stress,
						  const int stressSize,
						  const double* properties,
						  const int numProperties,
						  const double* totalStrain,
						  const int strainSize,
						  const bool computeStateVars)
{ // _calcStressElastic
  assert(0 != stress);
  assert(_MaxwellIsotropic3D::tensorSize == stressSize);
  assert(0 != properties);
  assert(_totalPropsQuadPt == numProperties);
  assert(0 != totalStrain);
  assert(_MaxwellIsotropic3D::tensorSize == strainSize);

  const double density = properties[_MaxwellIsotropic3D::pidDensity];
  const double mu = properties[_MaxwellIsotropic3D::pidMu];
  const double lambda = properties[_MaxwellIsotropic3D::pidLambda];
  const double maxwelltime = properties[_MaxwellIsotropic3D::pidMaxwellTime];

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
    // std::cout << "  " << properties[_MaxwellIsotropic3D::pidStrainT][iComp]
	    // << "   " << properties[_MaxwellIsotropic3D::pidVisStrain][iComp]
	    // << std::endl;

  const double meanStrainT = (properties[_MaxwellIsotropic3D::pidStrainT+0] +
			      properties[_MaxwellIsotropic3D::pidStrainT+1] +
			      properties[_MaxwellIsotropic3D::pidStrainT+2])/3.0;
  
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
    devStrainT = properties[_MaxwellIsotropic3D::pidStrainT+iComp] -
      diag[iComp]*meanStrainT;
    visStrain = expFac*properties[_MaxwellIsotropic3D::pidVisStrain+iComp] +
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
// Compute derivative of elasticity matrix at location from properties.
void
pylith::materials::MaxwellIsotropic3D::_calcElasticConstsElastic(
						  double* const elasticConsts,
						  const int numElasticConsts,
						  const double* properties,
						  const int numProperties,
						  const double* totalStrain,
						  const int strainSize)
{ // _calcElasticConstsElastic
  assert(0 != elasticConsts);
  assert(_MaxwellIsotropic3D::numElasticConsts == numElasticConsts);
  assert(0 != properties);
  assert(_totalPropsQuadPt == numProperties);
  assert(0 != totalStrain);
  assert(_MaxwellIsotropic3D::tensorSize == strainSize);
 
  const double density = properties[_MaxwellIsotropic3D::pidDensity];
  const double mu = properties[_MaxwellIsotropic3D::pidMu];
  const double lambda = properties[_MaxwellIsotropic3D::pidLambda];
  const double maxwelltime = properties[_MaxwellIsotropic3D::pidMaxwellTime];

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
// Compute derivative of elasticity matrix at location from properties
// as an elastic material.
void
pylith::materials::MaxwellIsotropic3D::_calcElasticConstsViscoelastic(
						  double* const elasticConsts,
						  const int numElasticConsts,
						  const double* properties,
						  const int numProperties,
						  const double* totalStrain,
						  const int strainSize)
{ // _calcElasticConstsViscoelastic
  assert(0 != elasticConsts);
  assert(_MaxwellIsotropic3D::numElasticConsts == numElasticConsts);
  assert(0 != properties);
  assert(_totalPropsQuadPt == numProperties);
  assert(0 != totalStrain);
  assert(_MaxwellIsotropic3D::tensorSize == strainSize);
 
  const double density = properties[_MaxwellIsotropic3D::pidDensity];
  const double mu = properties[_MaxwellIsotropic3D::pidMu];
  const double lambda = properties[_MaxwellIsotropic3D::pidLambda];
  const double maxwelltime = properties[_MaxwellIsotropic3D::pidMaxwellTime];

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
pylith::materials::MaxwellIsotropic3D::_updatePropertiesElastic(
						 double* const properties,
						 const int numProperties,
						 const double* totalStrain,
						 const int strainSize)
{ // _updatePropertiesElastic
  assert(0 != properties);
  assert(_totalPropsQuadPt == numProperties);
  assert(0 != totalStrain);
  assert(_MaxwellIsotropic3D::tensorSize == strainSize);

  const double maxwelltime = properties[_MaxwellIsotropic3D::pidMaxwellTime];

  const double e11 = totalStrain[0];
  const double e22 = totalStrain[1];
  const double e33 = totalStrain[2];

  const double traceStrainTpdt = e11 + e22 + e33;
  const double meanStrainTpdt = traceStrainTpdt/3.0;

  const double diag[] = { 1.0, 1.0, 1.0, 0.0, 0.0, 0.0 };

  // Temporary to get stresses.
  // double_array stress(6);
  // _calcStressElastic(&stress, (*properties), totalStrain);

  for (int iComp=0; iComp < _MaxwellIsotropic3D::tensorSize; ++iComp) {
    properties[_MaxwellIsotropic3D::pidStrainT+iComp] = totalStrain[iComp];
    properties[_MaxwellIsotropic3D::pidVisStrain+iComp] =
      totalStrain[iComp] - diag[iComp]*meanStrainTpdt;
  } // for
  PetscLogFlopsNoCheck(5 * _MaxwellIsotropic3D::tensorSize);
//   std::cout << std::endl;
//   std::cout << " updatePropertiesElastic: "<< std::endl;
//   std::cout << " StrainT  VisStrain  Stress: " << std::endl;
//   for (int iComp=0; iComp < _MaxwellIsotropic3D::tensorSize; ++iComp)
//     std::cout << "  " << properties[_MaxwellIsotropic3D::pidStrainT+iComp]
// 	    << "   " << properties[_MaxwellIsotropic3D::pidVisStrain+iComp]
// 	    << "   " << stress[iComp]
// 	    << std::endl;

  _needNewJacobian = true;
} // _updatePropertiesElastic

// ----------------------------------------------------------------------
// Update state variables.
void
pylith::materials::MaxwellIsotropic3D::_updatePropertiesViscoelastic(
						 double* const properties,
						 const int numProperties,
						 const double* totalStrain,
						 const int strainSize)
{ // _updatePropertiesViscoelastic
  assert(0 != properties);
  assert(_totalPropsQuadPt == numProperties);
  assert(0 != totalStrain);
  assert(_MaxwellIsotropic3D::tensorSize == strainSize);

  const double maxwelltime = properties[_MaxwellIsotropic3D::pidMaxwellTime];

  const double e11 = totalStrain[0];
  const double e22 = totalStrain[1];
  const double e33 = totalStrain[2];

  const double meanStrainTpdt = (e11 + e22 + e33) / 3.0;

  const double diag[] = { 1.0, 1.0, 1.0, 0.0, 0.0, 0.0 };

  // Temporary to get stresses.
  double stress[6];
  const int stressSize = 6;
  _calcStressViscoelastic(stress, stressSize, 
			  properties, numProperties,
			  totalStrain, strainSize, true);

  const double meanStrainT = 
    (properties[_MaxwellIsotropic3D::pidStrainT+0] +
     properties[_MaxwellIsotropic3D::pidStrainT+1] +
     properties[_MaxwellIsotropic3D::pidStrainT+2]) / 3.0;
  
  PetscLogFlopsNoCheck(6);

  double dq = ViscoelasticMaxwell::computeVisStrain(_dt, maxwelltime);

  const double expFac = exp(-_dt/maxwelltime);
  double devStrainTpdt = 0.0;
  double devStrainT = 0.0;
  double visStrain = 0.0;
  PetscLogFlopsNoCheck(3);
  for (int iComp=0; iComp < _MaxwellIsotropic3D::tensorSize; ++iComp) {
    devStrainTpdt = totalStrain[iComp] - diag[iComp]*meanStrainTpdt;
    devStrainT = properties[_MaxwellIsotropic3D::pidStrainT+iComp] -
      diag[iComp] * meanStrainT;
    visStrain = expFac * 
      properties[_MaxwellIsotropic3D::pidVisStrain+iComp] +
      dq * (devStrainTpdt - devStrainT);
    properties[_MaxwellIsotropic3D::pidVisStrain+iComp] = visStrain;
    properties[_MaxwellIsotropic3D::pidStrainT+iComp] = totalStrain[iComp];
  } // for
  PetscLogFlopsNoCheck(8 * _MaxwellIsotropic3D::tensorSize);

  _needNewJacobian = false;

//   std::cout << std::endl;
//   std::cout << " updatePropertiesViscoelastic: "<< std::endl;
//   std::cout << " StrainT  VisStrain  Stress: " << std::endl;
//   for (int iComp=0; iComp < _MaxwellIsotropic3D::tensorSize; ++iComp)
//     std::cout << "  " << properties[_MaxwellIsotropic3D::pidStrainT+iComp]
// 	    << "   " << properties[_MaxwellIsotropic3D::pidVisStrain+iComp]
// 	    << "   " << stress[iComp]
// 	    << std::endl;
} // _updatePropertiesViscoelastic


// End of file 
