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

#include "ElasticPlaneStrain.hh" // implementation of object methods

#include "petsc.h" // USES PetscLogFlopsNoCheck

#include <assert.h> // USES assert()
#include <sstream> // USES std::ostringstream
#include <stdexcept> // USES std::runtime_error

// ----------------------------------------------------------------------
namespace pylith {
  namespace materials {
    namespace _ElasticPlaneStrain {

      /// Number of entries in stress tensor.
      const int tensorSize = 3;

      /// Number of elastic constants (for general 3-D elastic material)
      const int numElasticConsts = 6;

      /// Number of physical properties.
      const int numProperties = 3;

      /// Physical properties.
      const Material::PropMetaData properties[] = {
	{ "density", 1, SCALAR_FIELD },
	{ "mu", 1, SCALAR_FIELD },
	{ "lambda", 1, SCALAR_FIELD },
      };
      /// Indices of physical properties
      const int pidDensity = 0;
      const int pidMu = pidDensity + 1;
      const int pidLambda = pidMu + 1;

      /// Values expected in spatial database
      const int numDBValues = 3;
      const char* namesDBValues[] = { "density", "vs", "vp" };
      
      /// Indices of database values
      const int didDensity = 0;
      const int didVs = 1;
      const int didVp = 2;

    } // _ElasticPlaneStrain
  } // materials
} // pylith


// ----------------------------------------------------------------------
// Default constructor.
pylith::materials::ElasticPlaneStrain::ElasticPlaneStrain(void) :
  ElasticMaterial(_ElasticPlaneStrain::tensorSize,
		  _ElasticPlaneStrain::numElasticConsts,
		  _ElasticPlaneStrain::namesDBValues,
		  _ElasticPlaneStrain::numDBValues,
		  _ElasticPlaneStrain::properties,
		  _ElasticPlaneStrain::numProperties)
{ // constructor
  _dimension = 2;
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::materials::ElasticPlaneStrain::~ElasticPlaneStrain(void)
{ // destructor
} // destructor

// ----------------------------------------------------------------------
// Compute parameters from values in spatial database.
void
pylith::materials::ElasticPlaneStrain::_dbToProperties(
				          double* const propValues,
                                          const double_array& dbValues) const
{ // _dbToProperties
  assert(0 != propValues);
  const int numDBValues = dbValues.size();
  assert(_ElasticPlaneStrain::numDBValues == numDBValues);

  const double density = dbValues[_ElasticPlaneStrain::didDensity];
  const double vs = dbValues[_ElasticPlaneStrain::didVs];
  const double vp = dbValues[_ElasticPlaneStrain::didVp];
 
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
  
  propValues[_ElasticPlaneStrain::pidDensity] = density;
  propValues[_ElasticPlaneStrain::pidMu] = mu;
  propValues[_ElasticPlaneStrain::pidLambda] = lambda;

  PetscLogFlopsNoCheck(6);
} // _dbToProperties

// ----------------------------------------------------------------------
// Compute density at location from properties.
void
pylith::materials::ElasticPlaneStrain::_calcDensity(
				  double* const density,
				  const double* properties,
				  const int numProperties)
{ // calcDensity
  assert(0 != density);
  assert(0 != properties);
  assert(_totalPropsQuadPt == numProperties);

  density[0] = properties[_ElasticPlaneStrain::pidDensity];
} // calcDensity

// ----------------------------------------------------------------------
// Compute stress tensor at location from properties.
void
pylith::materials::ElasticPlaneStrain::_calcStress(
				  double* const stress,
				  const int stressSize,
				  const double* properties,
				  const int numProperties,
				  const double* totalStrain,
				  const int strainSize,
				  const bool computeStateVars)
{ // _calcStress
  assert(0 != stress);
  assert(_ElasticPlaneStrain::tensorSize == stressSize);
  assert(0 != properties);
  assert(_totalPropsQuadPt == numProperties);
  assert(0 != totalStrain);
  assert(_ElasticPlaneStrain::tensorSize == strainSize);

  const double density = properties[_ElasticPlaneStrain::pidDensity];
  const double mu = properties[_ElasticPlaneStrain::pidMu];
  const double lambda = properties[_ElasticPlaneStrain::pidLambda];
  
  const double mu2 = 2.0 * mu;
  
  const double e11 = totalStrain[0];
  const double e22 = totalStrain[1];
  const double e12 = totalStrain[2];

  const double s12 = lambda * (e11 + e22);

  stress[0] = s12 + mu2*e11;
  stress[1] = s12 + mu2*e22;
  stress[2] = mu2 * e12;

  PetscLogFlopsNoCheck(8);
} // _calcStress

// ----------------------------------------------------------------------
// Compute density at location from properties.
void
pylith::materials::ElasticPlaneStrain::_calcElasticConsts(
					       double* const elasticConsts,
					       const int numElasticConsts,
					       const double* properties,
					       const int numProperties,
					       const double* totalStrain,
					       const int strainSize)
{ // calcElasticConsts
  assert(0 != elasticConsts);
  assert(_ElasticPlaneStrain::numElasticConsts == numElasticConsts);
  assert(0 != properties);
  assert(_totalPropsQuadPt == numProperties);
  assert(0 != totalStrain);
  assert(_ElasticPlaneStrain::tensorSize == strainSize);

  const double density = properties[_ElasticPlaneStrain::pidDensity];
  const double mu = properties[_ElasticPlaneStrain::pidMu];
  const double lambda = properties[_ElasticPlaneStrain::pidLambda];

  const double mu2 = 2.0 * mu;
  const double lambda2mu = lambda + mu2;
  
  elasticConsts[0] = lambda2mu; // C1111
  elasticConsts[1] = lambda; // C1122
  elasticConsts[2] = 0; // C1112
  elasticConsts[3] = lambda2mu; // C2222
  elasticConsts[4] = 0; // C2212
  elasticConsts[5] = mu2; // C1212

  PetscLogFlopsNoCheck(2);
} // calcElasticConsts


// End of file 
