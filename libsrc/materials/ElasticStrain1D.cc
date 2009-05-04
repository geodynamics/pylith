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

#include "ElasticStrain1D.hh" // implementation of object methods

#include "Metadata.hh" // USES Metadata

#include "pylith/utils/array.hh" // USES double_array
#include "pylith/utils/constdefs.h" // USES MAXDOUBLE

#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include "petsc.h" // USES PetscLogFlops

#include <cassert> // USES assert()
#include <sstream> // USES std::ostringstream
#include <stdexcept> // USES std::runtime_error

// ----------------------------------------------------------------------
namespace pylith {
  namespace materials {
    namespace _ElasticStrain1D {

      // Dimension of material.
      const int dimension = 1;

      /// Number of entries in stress tensor.
      const int tensorSize = 1;

      /// Number of entries in derivative of elasticity matrix.
      const int numElasticConsts = 1;

      /// Number of physical properties.
      const int numProperties = 3;

      // Physical properties.
      const Metadata::ParamDescription properties[] = {
	{ "density", 1, pylith::topology::FieldBase::SCALAR },
	{ "mu", 1, pylith::topology::FieldBase::SCALAR },
	{ "lambda", 1, pylith::topology::FieldBase::SCALAR },
      };

      // Values expected in spatial database
      const int numDBProperties = 3;
      const char* dbProperties[] = { "density", "vs", "vp" };      
      
    } // _ElasticStrain1D
  } // materials
} // pylith

// Indices of physical properties
const int pylith::materials::ElasticStrain1D::p_density = 0;

const int pylith::materials::ElasticStrain1D::p_mu = 
  pylith::materials::ElasticStrain1D::p_density + 1;

const int pylith::materials::ElasticStrain1D::p_lambda = 
  pylith::materials::ElasticStrain1D::p_mu + 1;

// Indices of database values (order must match dbProperties)
const int pylith::materials::ElasticStrain1D::db_density = 0;

const int pylith::materials::ElasticStrain1D::db_vs = 
  pylith::materials::ElasticStrain1D::db_density + 1;

const int pylith::materials::ElasticStrain1D::db_vp = 
  pylith::materials::ElasticStrain1D::db_vs + 1;

// ----------------------------------------------------------------------
// Default constructor.
pylith::materials::ElasticStrain1D::ElasticStrain1D(void) :
  ElasticMaterial(_ElasticStrain1D::dimension,
		  _ElasticStrain1D::tensorSize,
		  _ElasticStrain1D::numElasticConsts,
		  Metadata(_ElasticStrain1D::properties,
			   _ElasticStrain1D::numProperties,
			   _ElasticStrain1D::dbProperties,
			   _ElasticStrain1D::numDBProperties,
			   0, 0,
			   0, 0))
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::materials::ElasticStrain1D::~ElasticStrain1D(void)
{ // destructor
} // destructor

// ----------------------------------------------------------------------
// Compute parameters from values in spatial database.
void
pylith::materials::ElasticStrain1D::_dbToProperties(
					    double* const propValues,
					    const double_array& dbValues) const
{ // _dbToProperties
  assert(0 != propValues);
  const int numDBValues = dbValues.size();
  assert(_ElasticStrain1D::numDBProperties == numDBValues);

  const double density = dbValues[db_density];
  const double vs = dbValues[db_vs];
  const double vp = dbValues[db_vp];
 
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

  propValues[p_density] = density;
  propValues[p_mu] = mu;
  propValues[p_lambda] = lambda;

  PetscLogFlops(6);
} // _dbToProperties

// ----------------------------------------------------------------------
// Nondimensionalize properties.
void
pylith::materials::ElasticStrain1D::_nondimProperties(double* const values,
							 const int nvalues) const
{ // _nondimProperties
  assert(0 != _normalizer);
  assert(0 != values);
  assert(nvalues == _ElasticStrain1D::numProperties);

  const double densityScale = _normalizer->densityScale();
  const double pressureScale = _normalizer->pressureScale();

  values[p_density] = 
    _normalizer->nondimensionalize(values[p_density], densityScale);
  values[p_mu] = 
    _normalizer->nondimensionalize(values[p_mu], pressureScale);
  values[p_lambda] = 
    _normalizer->nondimensionalize(values[p_lambda], pressureScale);

  PetscLogFlops(3);
} // _nondimProperties

// ----------------------------------------------------------------------
// Dimensionalize properties.
void
pylith::materials::ElasticStrain1D::_dimProperties(double* const values,
						      const int nvalues) const
{ // _dimProperties
  assert(0 != _normalizer);
  assert(0 != values);
  assert(nvalues == _ElasticStrain1D::numProperties);

  const double densityScale = _normalizer->densityScale();
  const double pressureScale = _normalizer->pressureScale();

  values[p_density] = 
    _normalizer->dimensionalize(values[p_density], densityScale);
  values[p_mu] = 
    _normalizer->dimensionalize(values[p_mu], pressureScale);
  values[p_lambda] = 
    _normalizer->dimensionalize(values[p_lambda], pressureScale);

  PetscLogFlops(3);
} // _dimProperties

// ----------------------------------------------------------------------
// Compute density at location from properties.
void
pylith::materials::ElasticStrain1D::_calcDensity(double* const density,
						 const double* properties,
						 const int numProperties,
						 const double* stateVars,
						 const int numStateVars)
{ // _calcDensity
  assert(0 != density);
  assert(0 != properties);
  assert(_numPropsQuadPt == numProperties);
  assert(0 == numStateVars);

  density[0] = properties[p_density];
} // _calcDensity

// ----------------------------------------------------------------------
// Compute stress tensor at location from properties.
void
pylith::materials::ElasticStrain1D::_calcStress(double* const stress,
						const int stressSize,
						const double* properties,
						const int numProperties,
						const double* stateVars,
						const int numStateVars,
						const double* totalStrain,
						const int strainSize,
						const double* initialStress,
						const int initialStressSize,
						const double* initialStrain,
						const int initialStrainSize,
						const bool computeStateVars)
{ // _calcStress
  assert(0 != stress);
  assert(_ElasticStrain1D::tensorSize == stressSize);
  assert(0 != properties);
  assert(_numPropsQuadPt == numProperties);
  assert(0 == numStateVars);
  assert(0 != totalStrain);
  assert(_ElasticStrain1D::tensorSize == strainSize);
  assert(0 != initialStress);
  assert(_ElasticStrain1D::tensorSize == initialStressSize);
  assert(0 != initialStrain);
  assert(_ElasticStrain1D::tensorSize == initialStrainSize);

  const double density = properties[p_density];
  const double mu = properties[p_mu];
  const double lambda = properties[p_lambda];

  const double e11 = totalStrain[0] - initialStrain[0];
  stress[0] = (lambda + 2.0*mu) * e11 + initialStress[0];

  PetscLogFlops(5);
} // _calcStress

// ----------------------------------------------------------------------
// Compute derivative of elasticity matrix at location from properties.
void
pylith::materials::ElasticStrain1D::_calcElasticConsts(
					     double* const elasticConsts,
					     const int numElasticConsts,
					     const double* properties,
					     const int numProperties,
					     const double* stateVars,
					     const int numStateVars,
					     const double* totalStrain,
					     const int strainSize,
					     const double* initialStress,
					     const int initialStressSize,
					     const double* initialStrain,
					     const int initialStrainSize)
{ // _calcElasticConsts
  assert(0 != elasticConsts);
  assert(_ElasticStrain1D::numElasticConsts == numElasticConsts);
  assert(0 != properties);
  assert(_numPropsQuadPt == numProperties);
  assert(0 == numStateVars);
  assert(0 != totalStrain);
  assert(_ElasticStrain1D::tensorSize == strainSize);
  assert(0 != initialStress);
  assert(_ElasticStrain1D::tensorSize == initialStressSize);
  assert(0 != initialStrain);
  assert(_ElasticStrain1D::tensorSize == initialStrainSize);
 
  const double density = properties[p_density];
  const double mu = properties[p_mu];
  const double lambda = properties[p_lambda];

  elasticConsts[0] = lambda + 2.0*mu;

  PetscLogFlops(2);
} // _calcElasticConsts

// ----------------------------------------------------------------------
// Get stable time step for implicit time integration.
double
pylith::materials::ElasticStrain1D::_stableTimeStepImplicit(
					      const double* properties,
					      const int numProperties,
					      const double* stateVars,
					      const int numStateVars) const
{ // _stableTimeStepImplicit
  return pylith::PYLITH_MAXDOUBLE;
} // _stableTimeStepImplicit


// End of file 
