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

#include "ElasticStress1D.hh" // implementation of object methods

#include "petsc.h" // USES PetscLogFlopsNoCheck

#include <assert.h> // USES assert()

// ----------------------------------------------------------------------
namespace pylith {
  namespace materials {
    namespace _ElasticStress1D {

      /// Number of entries in stress tensor.
      const int tensorSize = 1;

      /// Number of entries in derivative of elasticity matrix.
      const int numElasticConsts = 1;

      /// Number of physical properties.
      const int numProperties = 3;

      /// Physical properties.
      const Material::PropMetaData properties[] = {
	{ "density", 1, OTHER_FIELD },
	{ "mu", 1, OTHER_FIELD },
	{ "lambda", 1, OTHER_FIELD },
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
    } // _ElasticStress1D
  } // materials
} // pylith

// ----------------------------------------------------------------------
// Default constructor.
pylith::materials::ElasticStress1D::ElasticStress1D(void) :
  ElasticMaterial(_ElasticStress1D::tensorSize,
		  _ElasticStress1D::numElasticConsts,
		  _ElasticStress1D::namesDBValues,
		  _ElasticStress1D::numDBValues,
		  _ElasticStress1D::properties,
		  _ElasticStress1D::numProperties)
{ // constructor
  _dimension = 1;
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::materials::ElasticStress1D::~ElasticStress1D(void)
{ // destructor
} // destructor

// ----------------------------------------------------------------------
// Compute parameters from values in spatial database.
void
pylith::materials::ElasticStress1D::_dbToProperties(
					      double* const propValues,
					      const double_array& dbValues) const
{ // _dbToProperties
  assert(0 != propValues);
  const int numDBValues = dbValues.size();
  assert(_ElasticStress1D::numDBValues == numDBValues);

  const double density = dbValues[_ElasticStress1D::didDensity];
  const double vs = dbValues[_ElasticStress1D::didVs];
  const double vp = dbValues[_ElasticStress1D::didVp];
 
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

  propValues[_ElasticStress1D::pidDensity] = density;
  propValues[_ElasticStress1D::pidMu] = mu;
  propValues[_ElasticStress1D::pidLambda] = lambda;

  PetscLogFlopsNoCheck(6);
} // _dbToProperties

// ----------------------------------------------------------------------
// Compute density at location from properties.
void
pylith::materials::ElasticStress1D::_calcDensity(double* const density,
						 const double* properties,
						 const int numProperties)
{ // _calcDensity
  assert(0 != density);
  assert(0 != properties);
  assert(_totalPropsQuadPt == numProperties);

  density[0] = properties[_ElasticStress1D::pidDensity];
} // _calcDensity

// ----------------------------------------------------------------------
// Compute stress tensor at location from properties.
void
pylith::materials::ElasticStress1D::_calcStress(double* const stress,
						const int stressSize,
						const double* properties,
						const int numProperties,
						const double* totalStrain,
						const int strainSize,
						const bool computeStateVars)
{ // _calcStress
  assert(0 != stress);
  assert(_ElasticStress1D::tensorSize == stressSize);
  assert(0 != properties);
  assert(_totalPropsQuadPt == numProperties);
  assert(0 != totalStrain);
  assert(_ElasticStress1D::tensorSize == strainSize);

  const double density = properties[_ElasticStress1D::pidDensity];
  const double mu = properties[_ElasticStress1D::pidMu];
  const double lambda = properties[_ElasticStress1D::pidLambda];

  const double e11 = totalStrain[0];
  stress[0] = mu * (3.0*lambda+2.0*mu) / (lambda + mu) * e11;

  PetscLogFlopsNoCheck(7);
} // _calcStress

// ----------------------------------------------------------------------
// Compute derivative of elasticity matrix at location from properties.
void
pylith::materials::ElasticStress1D::_calcElasticConsts(
				  double* const elasticConsts,
				  const int numElasticConsts,
				  const double* properties,
				  const int numProperties,
				  const double* totalStrain,
				  const int strainSize)
{ // _calcElasticConsts
  assert(0 != elasticConsts);
  assert(_ElasticStress1D::numElasticConsts == numElasticConsts);
  assert(0 != properties);
  assert(_totalPropsQuadPt == numProperties);
  assert(0 != totalStrain);
  assert(_ElasticStress1D::tensorSize == strainSize);
 
  const double density = properties[_ElasticStress1D::pidDensity];
  const double mu = properties[_ElasticStress1D::pidMu];
  const double lambda = properties[_ElasticStress1D::pidLambda];

  elasticConsts[0] = mu * (3.0*lambda+2.0*mu) / (lambda + mu);

  PetscLogFlopsNoCheck(6);
} // _calcElasticConsts


// End of file 
