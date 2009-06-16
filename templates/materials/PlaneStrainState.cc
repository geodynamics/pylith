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

#include "PlaneStrainState.hh" // implementation of object methods

#include "pylith/materials/Metadata.hh" // USES Metadata

#include "pylith/utils/array.hh" // USES double_array
#include "pylith/utils/constdefs.h" // USES MAXDOUBLE

#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

//#include "petsc.h" // USES PetscLogFlops

#include <cassert> // USES assert()
#include <sstream> // USES std::ostringstream
#include <stdexcept> // USES std::runtime_error

// ----------------------------------------------------------------------
// Create a local namespace to use for local constants and other
// information. This insulates all other classes from this information
// while preventing clashes with other local constants and data (as
// long as no other object use the same _PlaneStrainState namespace in
// the contrib::materials namespace.
namespace contrib {
  namespace materials {
    namespace _PlaneStrainState {

      // Dimension of material.
      const int dimension = 2;

      // Number of entries in stress tensor.
      const int tensorSize = 3;

      // Number of elastic constants (for general 3-D elastic material)
      const int numElasticConsts = 6;

      // These are the physical properties stored during the
      // simulation and need not coincide with the physical properties
      // provided by the user. We store Lame's constants and density,
      // but the user provides Vs, Vp, and density.

      // Number of physical properties.
      const int numProperties = 3;

      // Physical properties. 
      const pylith::materials::Metadata::ParamDescription properties[] = {
	{ "density", 1, pylith::topology::FieldBase::SCALAR },
	{ "mu", 1, pylith::topology::FieldBase::SCALAR },
	{ "lambda", 1, pylith::topology::FieldBase::SCALAR },
      };

      // Values expected in spatial database
      const int numDBProperties = 3;
      const char* dbProperties[] = { "density", "vs", "vp" };      
      
      // These are the state variables stored during the
      // simulation. Usually, we store only the time-dependent values
      // needed to compute the behavior at a given point in time. In
      // this example, however, for illustration purposes we store the
      // elastic strain (total strain) tensor and the stress tensor.

      /// Number of state variables. we store 
      const int numStateVars = 2;
      
      /// State variables.
      const pylith::materials::Metadata::ParamDescription stateVars[] = {
	{ "total_strain", tensorSize, pylith::topology::FieldBase::TENSOR },
	{ "stress", tensorSize, pylith::topology::FieldBase::TENSOR },
      };

      // We do not include a list of values expected to be in the
      // initial state variables spatial database because the initial
      // state variables are the initial stress and strain, which are
      // already provided as special spatial databases.

    } // _PlaneStrainState
  } // materials
} // contrib

// Indices of physical properties in the properties array.
const int contrib::materials::PlaneStrainState::p_density = 0;

const int contrib::materials::PlaneStrainState::p_mu = 
  contrib::materials::PlaneStrainState::p_density + 1;

const int contrib::materials::PlaneStrainState::p_lambda = 
  contrib::materials::PlaneStrainState::p_mu + 1;

// Indices of database values (order must match dbProperties)
const int contrib::materials::PlaneStrainState::db_density = 0;

const int contrib::materials::PlaneStrainState::db_vs = 
  contrib::materials::PlaneStrainState::db_density + 1;

const int contrib::materials::PlaneStrainState::db_vp = 
  contrib::materials::PlaneStrainState::db_vs + 1;

// Indices of state variables in the state variables array.
const int contrib::materials::PlaneStrainState::s_totalStrain = 0;

const int contrib::materials::PlaneStrainState::s_stress = 
  contrib::materials::PlaneStrainState::s_totalStrain + 
  contrib::materials::_PlaneStrainState::tensorSize;

// ----------------------------------------------------------------------
// Default constructor.
contrib::materials::PlaneStrainState::PlaneStrainState(void) :
  pylith::materials::ElasticMaterial(_PlaneStrainState::dimension,
				     _PlaneStrainState::tensorSize,
				     _PlaneStrainState::numElasticConsts,
     pylith::materials::Metadata(_PlaneStrainState::properties,
				 _PlaneStrainState::numProperties,
				 _PlaneStrainState::dbProperties,
				 _PlaneStrainState::numDBProperties,
				 _PlaneStrainState::stateVars,
				 _PlaneStrainState::numStateVars,
				 0, 0))
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor.
contrib::materials::PlaneStrainState::~PlaneStrainState(void)
{ // destructor
} // destructor

// ----------------------------------------------------------------------
// Compute parameters from values in spatial database.
void
contrib::materials::PlaneStrainState::_dbToProperties(
				double* const propValues,
                                const pylith::double_array& dbValues) const
{ // _dbToProperties
  // Check consistency of arguments
  assert(0 != propValues);
  const int numDBValues = dbValues.size();
  assert(_PlaneStrainState::numDBProperties == numDBValues);

  // Extract values from array using our defined indices.
  const double density = dbValues[db_density];
  const double vs = dbValues[db_vs];
  const double vp = dbValues[db_vp];
 
  // Check for reasonable values. If user supplied unreasonable values
  // throw an exception.
  if (density <= 0.0 || vs <= 0.0 || vp <= 0.0) {
    std::ostringstream msg;
    msg << "Spatial database returned nonpositive value for physical "
	<< "properties.\n"
	<< "density: " << density << "\n"
	<< "vp: " << vp << "\n"
	<< "vs: " << vs << "\n";
    throw std::runtime_error(msg.str());
  } // if

  // Compute physical properties that we store from the user-supplied
  // physical properties.
  const double mu = density * vs*vs;
  const double lambda = density * vp*vp - 2.0*mu;

  // Check for reasonable values. If values are unreasonable throw an
  // exception.
  if (lambda <= 0.0) {
    std::ostringstream msg;
    msg << "Attempted to set Lame's constant lambda to nonpositive value.\n"
	<< "density: " << density << "\n"
	<< "vp: " << vp << "\n"
	<< "vs: " << vs << "\n";
    throw std::runtime_error(msg.str());
  } // if

  // Store computed physical properties in the properties array.
  propValues[p_density] = density;
  propValues[p_mu] = mu;
  propValues[p_lambda] = lambda;

  //  PetscLogFlops(6); // Add to logging of number of floating point operations
} // _dbToProperties

// ----------------------------------------------------------------------
// Nondimensionalize properties.
void
contrib::materials::PlaneStrainState::_nondimProperties(double* const values,
							 const int nvalues) const
{ // _nondimProperties
  // Check consistency of arguments.
  assert(0 != _normalizer);
  assert(0 != values);
  assert(nvalues == _PlaneStrainState::numProperties);

  // Get scales needed to nondimensional parameters from the
  // Nondimensional object.
  const double densityScale = _normalizer->densityScale();
  const double pressureScale = _normalizer->pressureScale();

  // Use the Nondimensional::nondimensionalize() function to
  // nondimensionalize the quantities using the appropriate scale.
  values[p_density] = 
    _normalizer->nondimensionalize(values[p_density], densityScale);
  values[p_mu] = 
    _normalizer->nondimensionalize(values[p_mu], pressureScale);
  values[p_lambda] = 
    _normalizer->nondimensionalize(values[p_lambda], pressureScale);

  //  PetscLogFlops(3); // Add to logging of number of floating point operations
} // _nondimProperties

// ----------------------------------------------------------------------
// Dimensionalize properties.
void
contrib::materials::PlaneStrainState::_dimProperties(double* const values,
						      const int nvalues) const
{ // _dimProperties
  // Check consistency of arguments
  assert(0 != _normalizer);
  assert(0 != values);
  assert(nvalues == _PlaneStrainState::numProperties);

  // Get scales needed to dimensional parameters from the
  // Nondimensional object.
  const double densityScale = _normalizer->densityScale();
  const double pressureScale = _normalizer->pressureScale();

  // Use the Nondimensional::dimensionalize() function to
  // dimensionalize the quantities using the appropriate scale.
  values[p_density] = 
    _normalizer->dimensionalize(values[p_density], densityScale);
  values[p_mu] = 
    _normalizer->dimensionalize(values[p_mu], pressureScale);
  values[p_lambda] = 
    _normalizer->dimensionalize(values[p_lambda], pressureScale);

  //  PetscLogFlops(3); // Add to logging of number of floating point operations
} // _dimProperties

// ----------------------------------------------------------------------
// Compute density at location from properties.
void
contrib::materials::PlaneStrainState::_calcDensity(double* const density,
						    const double* properties,
						    const int numProperties,
						    const double* stateVars,
						    const int numStateVars)
{ // calcDensity
  // Check consistency of arguments.
  assert(0 != density);
  assert(0 != properties);
  assert(_numPropsQuadPt == numProperties);
  assert(0 != stateVars);
  assert(_numVarsQuadPt == numStateVars);

  // Set density using physical properties (trivial since one our
  // physical properties is density).
  density[0] = properties[p_density];
} // calcDensity

// ----------------------------------------------------------------------
// Compute stress tensor at location from properties.
void
contrib::materials::PlaneStrainState::_calcStress(double* const stress,
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
  // Check consistency of arguments.
  assert(0 != stress);
  assert(_PlaneStrainState::tensorSize == stressSize);
  assert(0 != properties);
  assert(_numPropsQuadPt == numProperties);
  assert(0 != stateVars);
  assert(_numVarsQuadPt == numStateVars);
  assert(0 != totalStrain);
  assert(_PlaneStrainState::tensorSize == strainSize);
  assert(0 != initialStress);
  assert(_PlaneStrainState::tensorSize == initialStressSize);
  assert(0 != initialStrain);
  assert(_PlaneStrainState::tensorSize == initialStrainSize);

  // Extract the material properties from the properties array.
  const double density = properties[p_density];
  const double mu = properties[p_mu];
  const double lambda = properties[p_lambda];

  // Compute a convenient constant.
  const double mu2 = 2.0*mu;

  // Compute the current strains accounting for the initial strain.
  const double e11 = totalStrain[0] - initialStrain[0];
  const double e22 = totalStrain[1] - initialStrain[1];
  const double e12 = totalStrain[2] - initialStrain[2];

  // Compute another convenient constant.
  const double s12 = lambda * (e11 + e22);

  // Compute the stresses and store them in the stress array.
  stress[0] = s12 + mu2*e11 + initialStress[0];
  stress[1] = s12 + mu2*e22 + initialStress[1];
  stress[2] = mu2 * e12 + initialStress[2];

  //  PetscLogFlops(14); // Add to logging of number of floating point operations
} // _calcStress

// ----------------------------------------------------------------------
// Compute elastic constants at location from properties.
void
contrib::materials::PlaneStrainState::_calcElasticConsts(
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
{ // calcElasticConsts
  // Check consistency of arguments.
  assert(0 != elasticConsts);
  assert(_PlaneStrainState::numElasticConsts == numElasticConsts);
  assert(0 != properties);
  assert(_numPropsQuadPt == numProperties);
  assert(0 != stateVars);
  assert(_numVarsQuadPt == numStateVars);
  assert(0 != totalStrain);
  assert(_PlaneStrainState::tensorSize == strainSize);
  assert(0 != initialStress);
  assert(_PlaneStrainState::tensorSize == initialStressSize);
  assert(0 != initialStrain);
  assert(_PlaneStrainState::tensorSize == initialStrainSize);
 
  // Extract the material properties from the properties array.
  const double density = properties[p_density];
  const double mu = properties[p_mu];
  const double lambda = properties[p_lambda];

  // Compute a couple convenient constants.
  const double mu2 = 2.0 * mu;
  const double lambda2mu = lambda + mu2;
   
  // Compute the elastic constants and store them in the elastic
  // constants array.
  elasticConsts[0] = lambda2mu; // C1111
  elasticConsts[1] = lambda; // C1122
  elasticConsts[2] = 0; // C1112
  elasticConsts[3] = lambda2mu; // C2222
  elasticConsts[4] = 0; // C2212
  elasticConsts[5] = mu2; // C1212

  //  PetscLogFlops(2); // Add to logging of number of floating point operations
} // calcElasticConsts

// ----------------------------------------------------------------------
// Get stable time step for implicit time integration.
double
contrib::materials::PlaneStrainState::stableTimeStepImplicit(
				    const pylith::topology::Mesh& mesh) {
  // Override the ElasticMaterial::stableTimeStepImplicit() function
  // (which calls _stableTimeStepImplicit() for each quadrature point
  // ) with an optimized calculation of the stable time step. This is
  // possible because the stable time step for an elastic material is
  // infinite, so we can simply return a very large number.
  return pylith::PYLITH_MAXDOUBLE;
}

// ----------------------------------------------------------------------
// Get stable time step for implicit time integration.
double
contrib::materials::PlaneStrainState::_stableTimeStepImplicit(
				     const double* properties,
				     const int numProperties,
				     const double* stateVars,
				     const int numStateVars) const
{ // _stableTimeStepImplicit Return the stable time step for this
  // material given its current state. This function will never be
  // called because we provide the stableTimeStepImplicit() function,
  // but we implement this function to satisfy the requirements of an
  // interface for an elastic material (which is defined by the
  // Material and ElasticMaterial objects).
  return pylith::PYLITH_MAXDOUBLE;
} // _stableTimeStepImplicit

// ----------------------------------------------------------------------
// Update state variables.
void
contrib::materials::PlaneStrainState::_updateStateVars(
					    double* const stateVars,
					    const int numStateVars,
					    const double* properties,
					    const int numProperties,
					    const double* totalStrain,
					    const int strainSize,
					    const double* initialStress,
					    const int initialStressSize,
					    const double* initialStrain,
					    const int initialStrainSize)
{ // _updateStateVars
  // Check consistency of arguments.
  assert(0 != stateVars);
  assert(_numVarsQuadPt == numStateVars);
  assert(0 != properties);
  assert(_numPropsQuadPt == numProperties);
  assert(0 != totalStrain);
  assert(_PlaneStrainState::tensorSize == strainSize);
  assert(0 != initialStress);
  assert(_PlaneStrainState::tensorSize == initialStressSize);
  assert(0 != initialStrain);
  assert(_PlaneStrainState::tensorSize == initialStrainSize);

  // Store the tensor size as a local value.
  const int tensorSize = _tensorSize;

  // Store the total strain in the state variable array.
  for (int iComp=0; iComp < tensorSize; ++iComp)
    stateVars[s_totalStrain+iComp] = totalStrain[iComp];

  // Use the _calcStress function to compute the stresses.
  const bool computeStateVars = false; // We are computing the state vars here!
  _calcStress(&stateVars[s_stress], tensorSize,
	      properties, numProperties,
	      stateVars, numStateVars,
	      totalStrain, strainSize,
	      initialStress, initialStressSize,
	      initialStrain, initialStrainSize,
	      computeStateVars);
} // _updateStateVars


// End of file 
