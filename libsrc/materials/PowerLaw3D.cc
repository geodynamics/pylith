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

#include "PowerLaw3D.hh" // implementation of object methods

#include "pylith/utils/array.hh" // USES double_array

#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include "petsc.h" // USES PetscLogFlops

#include <cassert> // USES assert()
#include <cstring> // USES memcpy()
#include <sstream> // USES std::ostringstream
#include <stdexcept> // USES std::runtime_error

// ----------------------------------------------------------------------
namespace pylith {
  namespace materials {
    namespace _PowerLaw3D{

      /// Dimension of material.
      const int dimension = 3;

      /// Number of entries in stress/strain tensors.
      const int tensorSize = 6;

      /// Number of entries in derivative of elasticity matrix.
      const int numElasticConsts = 21;

      /// Number of physical properties.
      const int numProperties = 5;

      /// Physical properties.
      const MetaData::ParamDescription properties[] = {
	{ "density", 1, pylith::topology::FieldBase::SCALAR },
	{ "mu", 1, pylith::topology::FieldBase::SCALAR },
	{ "lambda", 1, pylith::topology::FieldBase::SCALAR },
	{ "viscosity_coeff", 1, pylith::topology::FieldBase::SCALAR },
	{ "power_law_exponent", 1, pylith::topology::FieldBase::SCALAR }
      };

      // Values expected in properties spatial database
      const int numDBProperties = 5;
      const char* dbProperties[] = {"density", "vs", "vp" ,
				    "viscosity_coeff",
				    "power_law_exponent"};

      /// Number of state variables.
      const int numStateVars = 2;

      /// State variables.
      const Metadata::ParamDescription stateVars[] = {
	{ "viscous_strain", tensorSize, pylith::topology::FieldBase::TENSOR },
	{ "stress", tensorSize, pylith::topology::FieldBase::TENSOR }
      };

      // Values expected in state variables spatial database.
      const int numDBStateVars = 12;
      const char* dbStateVars[] = { "viscous-strain-xx",
				    "viscous-strain-yy",
				    "viscous-strain-zz",
				    "viscous-strain-xy",
				    "viscous-strain-yz",
				    "viscous-strain-xz",
				    "stress-xx",
				    "stress-yy",
				    "stress-zz",
				    "stress-xy",
				    "stress-yz",
				    "stress-xz"
      };

    } // _PowerLaw3D
  } // materials
} // pylith

// Indices of physical properties.
const int pylith::materials::PowerLaw3D::p_density = 0;

const int pylith::materials::PowerLaw3D::p_mu = 
  pylith::materials::PowerLaw3D::p_density + 1;

const int pylith::materials::PowerLaw3D::p_lambda = 
  pylith::materials::PowerLaw3D::p_mu + 1;

const int pylith::materials::PowerLaw3D::p_viscosityCoeff = 
  pylith::materials::PowerLaw3D::p_lambda + 1;

const int pylith::materials::PowerLaw3D::p_powerLawExponent = 
  pylith::materials::PowerLaw3D::p_viscosityCoeff + 1;

// Indices of property database values (order must match dbProperties).
const int pylith::materials::PowerLaw3D::db_density = 0;

const int pylith::materials::PowerLaw3D::db_vs = 
  pylith::materials::PowerLaw3D::db_density + 1;

const int pylith::materials::PowerLaw3D::db_vp = 
  pylith::materials::PowerLaw3D::db_vs + 1;

const int pylith::materials::PowerLaw3D::db_viscosityCoeff = 
  pylith::materials::PowerLaw3D::db_vp + 1;

const int pylith::materials::PowerLaw3D::db_powerLawExponent = 
  pylith::materials::PowerLaw3D::db_viscosityCoeff + 1;

// Indices of state variables.
const int pylith::materials::PowerLaw3D::s_viscousStrain = 0;

const int pylith::materials::PowerLaw3D::s_stress = 
  pylith::materials::PowerLaw3D::s_viscousStrain + 
  pylith::materials::PowerLaw3D::tensorSize;

// Indices of state variable database values (order must match dbStateVars).
const int pylith::materials::PowerLaw3D::db_viscousStrain = 0;

const int pylith::materials::PowerLaw3D::db_stress = 
  pylith::materials::PowerLaw3D::db_viscousStrain +
  pylith::materials::PowerLaw3D::tensorSize;

// ----------------------------------------------------------------------
// Default constructor.
pylith::materials::PowerLaw3D::PowerLaw3D(void) :
  ElasticMaterial(_PowerLaw3D::dimension,
		  _PowerLaw3D::tensorSize,
		  _PowerLaw3D::numElasticConsts,
		  Metadata(_PowerLaw3D::properties,
			   _PowerLaw3D::numProperties,
			   _PowerLaw3D::dbProperties,
			   _PowerLaw3D::numDBProperties,
			   _PowerLaw3D::stateVars,
			   _PowerLaw3D::numStateVars,
			   _PowerLaw3D::dbStateVars,
			   _PowerLaw3D::numDBStateVars)),
  _calcElasticConstsFn(0),
  _calcStressFn(0),
  _updateStateVarsFn(0)
{ // constructor
  useElasticBehavior(true);
  _viscousStrain.resize(_tensorSize);
  _stress.resize(_tensorSize);
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::materials::PowerLaw3D::~PowerLaw3D(void)
{ // destructor
} // destructor

// ----------------------------------------------------------------------
// Set whether elastic or inelastic constitutive relations are used.
void
pylith::materials::PowerLaw3D::useElasticBehavior(const bool flag)
{ // useElasticBehavior
  if (flag) {
    _calcStressFn = 
      &pylith::materials::PowerLaw3D::_calcStressElastic;
    _calcElasticConstsFn = 
      &pylith::materials::PowerLaw3D::_calcElasticConstsElastic;
    _updateStateVarsFn = 
      &pylith::materials::PowerLaw3D::_updateStateVarsElastic;

  } else {
    _calcStressFn = 
      &pylith::materials::PowerLaw3D::_calcStressViscoelastic;
    _calcElasticConstsFn = 
      &pylith::materials::PowerLaw3D::_calcElasticConstsViscoelastic;
    _updateStateVarsFn = 
      &pylith::materials::PowerLaw3D::_updateStateVarsViscoelastic;
  } // if/else
} // useElasticBehavior

// ----------------------------------------------------------------------
// Compute properties from values in spatial database.
void
pylith::materials::PowerLaw3D::_dbToProperties(
					    double* const propValues,
					    const double_array& dbValues) const
{ // _dbToProperties
  assert(0 != propValues);
  const int numDBValues = dbValues.size();
  assert(_PowerLaw3D::numDBValues == numDBValues);

  const double density = dbValues[db_density];
  const double vs = dbValues[db_vs];
  const double vp = dbValues[db_vp];
  const double viscosityCoeff = dbValues[db_viscosityCoeff];
  const double powerLawExponent = dbValues[db_powerLawExponent];
 
  if (density <= 0.0 || vs <= 0.0 || vp <= 0.0 || viscosityCoeff <= 0.0
      || powerLawExponent < 1.0) {
    std::ostringstream msg;
    msg << "Spatial database returned illegal value for physical "
	<< "properties.\n"
	<< "density: " << density << "\n"
	<< "vp: " << vp << "\n"
	<< "vs: " << vs << "\n"
	<< "viscosityCoeff: " << viscosityCoeff << "\n"
	<< "powerLawExponent: " << powerLawExponent << "\n";
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

  propValues[p_density] = density;
  propValues[p_mu] = mu;
  propValues[p_lambda] = lambda;
  propValues[p_viscosityCoeff] = viscosityCoeff;
  propValues[p_powerLawExponent] = powerLawExponent;

  PetscLogFlops(6);
} // _dbToProperties

// ----------------------------------------------------------------------
// Nondimensionalize properties.
void
pylith::materials::PowerLaw3D::_nondimProperties(double* const values,
							 const int nvalues) const
{ // _nondimProperties
  assert(0 != _normalizer);
  assert(0 != values);
  assert(nvalues == _numPropsQuadPt);

  const double densityScale = _normalizer->densityScale();
  const double pressureScale = _normalizer->pressureScale();
  const double timeScale = _normalizer->timeScale();
  // **** NOTE:  Make sure scaling is correct for viscosity coefficient.
  const double powerLawExponent = values[p_powerLawExponent];
  const double viscosityCoeffScale =
    (pressureScale^(1.0/powerLawExponent))/timeScale;
  values[p_density] = 
    _normalizer->nondimensionalize(values[p_density], densityScale);
  values[p_mu] = 
    _normalizer->nondimensionalize(values[p_mu], pressureScale);
  values[p_lambda] = 
    _normalizer->nondimensionalize(values[p_lambda], pressureScale);
  values[p_viscosityCoeff] = 
    _normalizer->nondimensionalize(values[p_viscosityCoeff],
				   viscosityCoeffScale);

  PetscLogFlops(7);
} // _nondimProperties

// ----------------------------------------------------------------------
// Dimensionalize properties.
void
pylith::materials::PowerLaw3D::_dimProperties(double* const values,
						      const int nvalues) const
{ // _dimProperties
  assert(0 != _normalizer);
  assert(0 != values);
  assert(nvalues == _totalPropsQuadPt);

  const double densityScale = _normalizer->densityScale();
  const double pressureScale = _normalizer->pressureScale();
  const double timeScale = _normalizer->timeScale();
  // **** NOTE:  Make sure scaling is correct for viscosity coefficient.
  const double powerLawExponent = values[p_powerLawExponent];
  const double viscosityCoeffScale =
    (pressureScale^(1.0/powerLawExponent))/timeScale;
  values[p_density] = 
    _normalizer->dimensionalize(values[p_density], densityScale);
  values[p_mu] = 
    _normalizer->dimensionalize(values[p_mu], pressureScale);
  values[p_lambda] = 
    _normalizer->dimensionalize(values[p_lambda], pressureScale);
  values[p_viscosityCoeff] = 
    _normalizer->dimensionalize(values[p_viscosityCoeff],
				viscosityCoeffScale);

  PetscLogFlops(7);
} // _dimProperties

// ----------------------------------------------------------------------
// Compute initial state variables from values in spatial database.
void
pylith::materials::PowerLaw3D::_dbToStateVars(
				double* const stateValues,
				const double_array& dbValues) const
{ // _dbToStateVars
  assert(0 != stateValues);
  const int numDBValues = dbValues.size();
  assert(_PowerLaw3D::numDBStateVars == numDBValues);

  const int totalSize = 2 * _tensorSize;
  assert(totalSize == _numVarsQuadPt);
  assert(totalSize == numDBValues);
  memcpy(stateValues[s_viscousStrain], &dbValues[db_viscousStrain],
	 _tensorSize*sizeof(double));
  memcpy(stateValues[s_stress], &dbValues[db_stress],
	 _tensorSize*sizeof(double));

  PetscLogFlops(0);
} // _dbToStateVars

// ----------------------------------------------------------------------
// Compute density at location from properties.
void
pylith::materials::PowerLaw3D::_calcDensity(double* const density,
					    const double* properties,
					    const int numProperties,
					    const double* stateVars,
					    const int numStateVars)
{ // _calcDensity
  assert(0 != density);
  assert(0 != properties);
  assert(_totalPropsQuadPt == numProperties);

  density[0] = properties[p_density];
} // _calcDensity

// ----------------------------------------------------------------------
// Get stable time step for implicit time integration.
double
pylith::materials::PowerLaw3D::_stableTimeStepImplicit(
				  const double* properties,
				  const int numProperties,
				  const double* stateVars,
				  const int numStateVars) const
{ // _stableTimeStepImplicit
  assert(0 != properties);
  assert(_numPropsQuadPt == numProperties);
  assert(0 != stateVars);
  assert(_numVarsQuadPt == numStateVars);

  memcpy(&stress[0], &stateVars[s_stress],
	   tensorSize * sizeof(double));
  const double meanStress = (stress[0] + stress[1] + stress[2])/3.0;
  const double devStress[] = {stress[0] - meanStress,
			      stress[1] - meanStress,
			      stress[2] - meanStress,
			      stress[3],
			      stress[4],
			      stress[5] };
  const double effStress = sqrt(0.5 * _scalarProduct(devStress, devStress));
  dtStable = 1.0;
  if (effStress != 0.0)
    dtStable = 0.1 * ((viscosityCoeff/effStressTpdt) ^ (powerLawExp - 1.0)) *
      (viscosityCoeff/mu);

  return dtStable;
} // _stableTimeStepImplicit

// ----------------------------------------------------------------------
// Compute stress tensor at location from properties as an elastic
// material.
void
pylith::materials::PowerLaw3D::_calcStressElastic(
				         double* const stress,
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
{ // _calcStressElastic
  assert(0 != stress);
  assert(_PowerLaw3D::tensorSize == stressSize);
  assert(0 != properties);
  assert(_numPropsQuadPt == numProperties);
  assert(0 != stateVars);
  assert(_numVarsQuadPt == numStateVars);
  assert(0 != totalStrain);
  assert(_PowerLaw3D::tensorSize == strainSize);
  assert(0 != initialStress);
  assert(_PowerLaw3D::tensorSize == initialStressSize);
  assert(0 != initialStrain);
  assert(_PowerLaw3D::tensorSize == initialStrainSize);

  const double mu = properties[p_mu];
  const double lambda = properties[p_lambda];
  const double mu2 = 2.0 * mu;

  const double e11 = totalStrain[0] - initialStrain[0];
  const double e22 = totalStrain[1] - initialStrain[1];
  const double e33 = totalStrain[2] - initialStrain[2];
  const double e12 = totalStrain[3] - initialStrain[3];
  const double e23 = totalStrain[4] - initialStrain[4];
  const double e13 = totalStrain[5] - initialStrain[5];
  
  const double traceStrainTpdt = e11 + e22 + e33;
  const double s123 = lambda * traceStrainTpdt;

  stress[0] = s123 + mu2*e11 + initialStress[0];
  stress[1] = s123 + mu2*e22 + initialStress[1];
  stress[2] = s123 + mu2*e33 + initialStress[2];
  stress[3] = mu2 * e12 + initialStress[3];
  stress[4] = mu2 * e23 + initialStress[4];
  stress[5] = mu2 * e13 + initialStress[5];

  PetscLogFlops(25);
} // _calcStressElastic

// ----------------------------------------------------------------------
// Compute stress tensor at location from properties as a viscoelastic
// material.
void
pylith::materials::PowerLaw3D::_calcStressViscoelastic(
					double* const stress,
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
{ // _calcStressViscoelastic
  assert(0 != stress);
  assert(_PowerLaw3D::tensorSize == stressSize);
  assert(0 != properties);
  assert(_numPropsQuadPt == numProperties);
  assert(0 != stateVars);
  assert(_numPropsQuadPt == numStateVars);
  assert(0 != totalStrain);
  assert(_PowerLaw3D::tensorSize == strainSize);
  assert(0 != initialStress);
  assert(_PowerLaw3D::tensorSize == initialStressSize);
  assert(0 != initialStrain);
  assert(_PowerLaw3D::tensorSize == initialStrainSize);

  const int tensorSize = _PowerLaw3D::tensorSize;
    
  // We need to do root-finding method if state variables are from previous
  // time step.
  if (computeStateVars) {

    const double mu = properties[p_mu];
    const double lambda = properties[p_lambda];
    const double viscosityCoeff = properties[p_viscosityCoeff];
    const double powerLawExp = properties[p_powerLawExponent];
    memcpy(&visStrainT[0], &stateVars[s_viscousStrain],
	   tensorSize * sizeof(double));
    memcpy(&stressT[0], &stateVars[s_stress],
	   tensorSize * sizeof(double));

    const double mu2 = 2.0 * mu;
    const double lamPlusMu = lambda + mu;
    const double bulkModulus = lambda + mu2/3.0;
    const double ae = 1.0/mu2;
    const double diag[] = { 1.0, 1.0, 1.0, 0.0, 0.0, 0.0 };

    // Need to figure out how time integration parameter alpha is going to be
    // specified.  It should probably be specified in the problem definition and
    // then used only by the material types that use it.  For now we are setting
    // it to 0.5, which should probably be the default value.
    const double alpha = 0.5;
    const double timeFac = _dt * (1.0 - alpha);

    // Initial stress values
    const double meanStressInitial = (initialStress[0] + initialStress[1] +
				      initialStress[2])/3.0;
    const double devStressInitial[] = { initialStress[0] - meanStressInitial,
					initialStress[1] - meanStressInitial,
					initialStress[2] - meanStressInitial,
					initialStress[3],
					initialStress[4],
					initialStress[5] };
    const double stressInvar2Initial = 0.5 *
      _scalarProduct(devStressInitial, devStressInitial);

    // Initial strain values
    const double meanStrainInitial = (initialStrain[0] + initialStrain[1] +
				      initialStrain[2])/3.0;

    // Values for current time step
    const double e11 = totalStrain[0];
    const double e22 = totalStrain[1];
    const double e33 = totalStrain[2];
    const double traceStrainTpdt = e11 + e22 + e33;
    const double meanStrainTpdt = traceStrainTpdt/3.0 - meanStrainInitial;
    const double meanStressTpdt = bulkModulus * traceStrainTpdt;

    // Note that I use the initial strain rather than the deviatoric initial
    // strain since otherwise the initial mean strain would get used twice.
    const double strainPPTpdt[] =
      { totalStrain[0] - meanStrainTpdt - visStrainT[0] - initialStrain[0],
	totalStrain[1] - meanStrainTpdt - visStrainT[1] - initialStrain[1],
	totalStrain[2] - meanStrainTpdt - visStrainT[2] - initialStrain[2],
	totalStrain[3] - visStrainT[3] - initialStrain[3],
	totalStrain[4] - visStrainT[4] - initialStrain[4],
	totalStrain[5] - visStrainT[5] - initialStrain[5] };
    const double strainPPInvar2Tpdt = 0.5 *
      _scalarProduct(strainPPTpdt, strainPPTpdt);

    // Values for previous time step
    const double meanStressT = (stressT[0] + stressT[1] + stressT[2])/3.0;
    const double devStressT[] = { stressT[0] - meanStressT,
				  stressT[1] - meanStressT,
				  stressT[2] - meanStressT,
				  stressT[3],
				  stressT[4],
				  stressT[5] };
    const double stressInvar2T = 0.5 * _scalarProduct(devStressT, devStressT);
    const double effStressT = sqrt(stressInvar2T);

    // Finish defining parameters needed for root-finding algorithm.
    const double b = strainInvar2Tpdt +
      ae * _scalarProduct(strainPPTpdt, devStressInitial) +
      ae * ae * stressInvar2Initial;
    const double c = (_scalarProduct(strainPPTpdt, devStressT) +
		      ae * _scalarProduct(devStressT, devStressInitial)) *
      timeFac;
    const double d = timeFac * effStressT;

    PetscLogFlops(55);
    // Put parameters into a vector and call root-finding algorithm.
    // This could also be a struct.
    const double effStressParams[] = {ae,
				      b,
				      c,
				      d,
				      alpha,
				      _dt,
				      effectiveStressT,
				      powerLawExp,
				      viscosityCoeff};
    // I think the PETSc root-finding procedure is too involved for what we want
    // here.  I would like the function to work something like:
    const double effStressInitialGuess = effStressT;
    double effStressTpdt =
      EffectiveStress::getEffStress(effStressInitialGuess,
				    effStressParams,
				    pylith::materials::PowerLaw3D::_effStressFunc,
				    pylith::materials::PowerLaw3D::_effStressFuncDFunc);

    // Compute stresses from effective stress.
    const double effStressTau = (1.0 - alpha) * effStressT +
      alpha * effStressTpdt;
    const double gammaTau = 0.5 *
      ((effStressTau/viscosityCoeff)^(powerLawExp - 1.0))/visscosityCoeff;
    const double factor1 = 1.0/(ae + alpha * _dt * gammaTau);
    const double factor2 = timeFac * gammaTau;
    double devStressTpdt = 0.0;

    for (int iComp=0; iComp < tensorSize; ++iComp) {
      devStressTpdt = factor1 *
	(strainPPTpdt[iComp] - factor2 * devStressT[iComp] +
	 ae * devStressInitial[iComp]);
      stress[iComp] = devStressTpdt + diag[iComp] *
	(meanStressTpdt + meanStressInitial);
    } // for
    PetscLogFlops(14 + 8 * tensorSize);

    // If state variables have already been updated, current stress is already
    // contained in stress.
  } else {
    memcpy(&stress[0], &stateVars[s_stress], tensorSize * sizeof(double));
  } // else

} // _calcStressViscoelastic

// ----------------------------------------------------------------------
// Effective stress function that computes effective stress function only
// (no derivative).
double
pylith::materials::PowerLaw3D::_effStressFunc(double effStressTpdt,
					      double *params)
{ // _effStressFunc
  double ae = params[0];
  double b = params[1];
  double c = params[2];
  double d = params[3];
  double alpha = params[4];
  double dt = params[5];
  double effStressT = params[6];
  double powerLawExp = params[7];
  double viscosityCoeff = params[8];
  double factor1 = 1.0-alpha;
  double effStressTau = factor1 * effStressT + alpha * effStressTpdt;
  double gammaTau = 0.5 * ((effStressTau/viscosityCoeff)^
			   (powerLawExp - 1.0))/viscosityCoeff;
  double a = ae + alpha * dt * gammaTau;
  double y = a * a * effStressTpdt * effStressTpdt - b +
    c * gammaTau - d * d * gammaTau * gammaTau;
  PetscLogFlops(21);
  return y;
} // _effStressFunc

// ----------------------------------------------------------------------
// Effective stress function that computes effective stress function
// derivative only (no function value).
double
pylith::materials::PowerLaw3D::_effStressDFunc(double effStressTpdt,
					       double *params)
{ // _effStressDFunc
  double ae = params[0];
  double c = params[2];
  double d = params[3];
  double alpha = params[4];
  double dt = params[5];
  double effStressT = params[6];
  double powerLawExp = params[7];
  double viscosityCoeff = params[8];
  double factor1 = 1.0-alpha;
  double effStressTau = factor1 * effStressT + alpha * effStressTpdt;
  double gammaTau = 0.5 * ((effStressTau/viscosityCoeff)^
			   (powerLawExp - 1.0))/viscosityCoeff;
  double a = ae + alpha * dt * gammaTau;
  double dGammaTau = 0.5 * alpha * (powerLawExp - 1.0) *
    ((effStressTau/viscosityCoeff) ^ (powerLawExp - 2.0))/
    (viscosityCoeff * viscosityCoeff);
  double dy = 2.0 * a * a * effStressTpdt + dGammaTau *
    (2.0 * a * alpha * dt * effStressTpdt * effStressTpdt +
     c - 2.0 * d * d * gammaTau);
  PetscLogFlops(36);
  return dy;
} // _effStressDFunc

// ----------------------------------------------------------------------
// Effective stress function that computes effective stress function
// and derivative.
void
pylith::materials::PowerLaw3D::_effStressFuncDFunc(double effStressTpdt,
						   double *params,
						   double *y,
						   double *dy)
{ // _effStressFunc
  double ae = params[0];
  double b = params[1];
  double c = params[2];
  double d = params[3];
  double alpha = params[4];
  double dt = params[5];
  double effStressT = params[6];
  double powerLawExp = params[7];
  double viscosityCoeff = params[8];
  double factor1 = 1.0-alpha;
  double effStressTau = factor1 * effStressT + alpha * effStressTpdt;
  double gammaTau = 0.5 * ((effStressTau/viscosityCoeff)^
			   (powerLawExp - 1.0))/viscosityCoeff;
  double dGammaTau = 0.5 * alpha * (powerLawExp - 1.0) *
    ((effStressTau/viscosityCoeff) ^ (powerLawExp - 2.0))/
    (viscosityCoeff * viscosityCoeff);
  double a = ae + alpha * dt * gammaTau;
  double *y = a * a * effStressTpdt * effStressTpdt - b +
    c * gammaTau - d * d * gammaTau * gammaTau;
  double *dy = 2.0 * a * a * effStressTpdt + dGammaTau *
    (2.0 * a * alpha * dt * effStressTpdt * effStressTpdt +
     c - 2.0 * d * d * gammaTau);
  PetscLogFlops(46);
} // _effStressFuncDFunc

// ----------------------------------------------------------------------
// Compute derivative of elasticity matrix at location from properties.
void
pylith::materials::PowerLaw3D::_calcElasticConstsElastic(
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
{ // _calcElasticConstsElastic
  assert(0 != elasticConsts);
  assert(_PowerLaw3D::numElasticConsts == numElasticConsts);
  assert(0 != properties);
  assert(_numPropsQuadPt == numProperties);
  assert(0 != stateVars);
  assert(_numVarsQuadPt == numStateVars);
  assert(0 != totalStrain);
  assert(_PowerLaw3D::tensorSize == strainSize);
  assert(0 != initialStress);
  assert(_PowerLaw3D::tensorSize == initialStressSize);
  assert(0 != initialStrain);
  assert(_PowerLaw3D::tensorSize == initialStrainSize);
 
  const double mu = properties[p_mu];
  const double lambda = properties[p_lambda];

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

  PetscLogFlops(2);
} // _calcElasticConstsElastic

// ----------------------------------------------------------------------
// Compute derivative of elasticity matrix at location from properties
// as a viscoelastic material. This version is to be used for the first
// iteration, before strains have been computed.
void
pylith::materials::PowerLaw3D::_calcElasticConstsViscoelasticInitial(
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
{ // _calcElasticConstsViscoelasticInitial
  assert(0 != elasticConsts);
  assert(_PowerLaw3D::numElasticConsts == numElasticConsts);
  assert(0 != properties);
  assert(_numPropsQuadPt == numProperties);
  assert(0 != stateVars);
  assert(_numPropsQuadPt == numStateVars);
  assert(0 != totalStrain);
  assert(_PowerLaw3D::tensorSize == strainSize);
  assert(0 != initialStress);
  assert(_PowerLaw3D::tensorSize == initialStressSize);
  assert(0 != initialStrain);
  assert(_PowerLaw3D::tensorSize == initialStrainSize);

  const int tensorSize = _PowerLaw3D::tensorSize;
  const double mu = properties[p_mu];
  const double lambda = properties[p_lambda];
  const double viscosityCoeff = properties[p_viscosityCoeff];
  const double powerLawExp = properties[p_powerLawExponent];
  memcpy(&stress[0], &stateVars[s_stress], tensorSize * sizeof(double));

  const double mu2 = 2.0 * mu;
  const double ae = 1.0/mu2;
  const double bulkModulus = lambda + mu2/3.0;

  const double meanStress = (stress[0] + stress[1] + stress[2])/3.0;
  const double devStress[] = {stress[0] - meanStress,
			      stress[1] - meanStress,
			      stress[2] - meanStress,
			      stress[3],
			      stress[4],
			      stress[5]};
  const double effStress = sqrt(0.5 * _scalarProduct(devStress, devStress));
  const double gamma = 0.5 *
    ((viscosityCoeff/effStress)^(powerLawExp - 1.0))/viscosityCoeff;
  const double visFac = 1.0/(3.0 * (ae + _dt * gamma));

  elasticConsts[ 0] = bulkModulus + 2.0 * visFac; // C1111
  elasticConsts[ 1] = bulkModulus - visFac; // C1122
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
  elasticConsts[15] = 3.0 * visFac; // C1212
  elasticConsts[16] = 0; // C1223
  elasticConsts[17] = 0; // C1213
  elasticConsts[18] = elasticConsts[15]; // C2323
  elasticConsts[19] = 0; // C2313
  elasticConsts[20] = elasticConsts[15]; // C1313

  PetscLogFlops(25);
} // _calcElasticConstsViscoelasticInitial

// ----------------------------------------------------------------------
// Compute derivative of elasticity matrix at location from properties
// as a viscoelastic material. This version is to be used after the first
// iteration, once strains have already been computed.
void
pylith::materials::PowerLaw3D::_calcElasticConstsViscoelastic(
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
{ // _calcElasticConstsViscoelastic
  assert(0 != elasticConsts);
  assert(_PowerLaw3D::numElasticConsts == numElasticConsts);
  assert(0 != properties);
  assert(_numPropsQuadPt == numProperties);
  assert(0 != stateVars);
  assert(_numPropsQuadPt == numStateVars);
  assert(0 != totalStrain);
  assert(_PowerLaw3D::tensorSize == strainSize);
  assert(0 != initialStress);
  assert(_PowerLaw3D::tensorSize == initialStressSize);
  assert(0 != initialStrain);
  assert(_PowerLaw3D::tensorSize == initialStrainSize);

  const double mu = properties[p_mu];
  const double lambda = properties[p_lambda];
  const double viscosityCoeff = properties[p_viscosityCoeff];
  const double powerLawExp = properties[p_powerLawExponent];
  memcpy(&visStrainT[0], &stateVars[s_viscousStrain],
	 tensorSize * sizeof(double));
  memcpy(&stressT[0], &stateVars[s_stress], tensorSize * sizeof(double));

  const double mu2 = 2.0 * mu;
  const double lamPlusMu = lambda + mu;
  const double bulkModulus = lambda + mu2/3.0;
  const double ae = 1.0/mu2;
  const double diag[] = { 1.0, 1.0, 1.0, 0.0, 0.0, 0.0 };
  
  // Need to figure out how time integration parameter alpha is going to be
  // specified.  It should probably be specified in the problem definition and
  // then used only by the material types that use it.  For now we are setting
  // it to 0.5, which should probably be the default value.
  const double alpha = 0.5;
  const double timeFac = _dt * (1.0 - alpha);
  
  // Values for current time step
  const double e11 = totalStrain[0];
  const double e22 = totalStrain[1];
  const double e33 = totalStrain[2];
  const double traceStrainTpdt = e11 + e22 + e33;
  const double meanStrainTpdt = traceStrainTpdt/3.0;
  const double meanStressTpdt = bulkModulus * traceStrainTpdt;

  // Initial stress values
  const double meanStressInitial = (_stressInitial[0] + _stressInitial[1] +
				    _stressInitial[2])/3.0;
  const double devStressInitial[] = { _stressInitial[0] - meanStressInitial,
				      _stressInitial[1] - meanStressInitial,
				      _stressInitial[2] - meanStressInitial,
				      _stressInitial[3],
				      _stressInitial[4],
				      _stressInitial[5] };
  const double stressInvar2Initial = 0.5 *
    _scalarProduct(devStressInitial, devStressInitial);
  
  // Values for current time step
  const double strainPPTpdt[] =
    { _totalStrain[0] - meanStrainTpdt - visStrainT[0],
      _totalStrain[1] - meanStrainTpdt - visStrainT[1],
      _totalStrain[2] - meanStrainTpdt - visStrainT[2],
      _totalStrain[3] - visStrainT[3],
      _totalStrain[4] - visStrainT[4],
      _totalStrain[5] - visStrainT[5] };
  const double strainPPInvar2Tpdt = 0.5 *
    _scalarProduct(strainPPTpdt, strainPPTpdt);
  
  // Values for previous time step
  const double meanStressT = (_stressT[0] + _stressT[1] + _stressT[2])/3.0;
  const double devStressT[] = { _stressT[0] - meanStressT,
				_stressT[1] - meanStressT,
				_stressT[2] - meanStressT,
				_stressT[3],
				_stressT[4],
				_stressT[5] };
  const double stressInvar2T = 0.5 * _scalarProduct(devStressT, devStressT);
  const double effStressT = sqrt(stressInvar2T);
  
  // Finish defining parameters needed for root-finding algorithm.
  const double b = strainInvar2Tpdt +
    ae * _scalarProduct(strainPPTpdt, devStressInitial) +
    ae * ae * stressInvar2Initial;
  const double c = (_scalarProduct(strainPPTpdt, devStressT) +
		    ae * _scalarProduct(devStressT, devStressInitial)) *
    timeFac;
  const double d = timeFac * effStressT;
  
  PetscLogFlops(45);
  // Put parameters into a vector and call root-finding algorithm.
  // This could also be a struct.
  const double effStressParams[] = {ae,
				    b,
				    c,
				    d,
				    alpha,
				    _dt,
				    effectiveStressT,
				    powerLawExp,
				    viscosityCoeff};
  // I think the PETSc root-finding procedure is too involved for what we want
  // here.  I would like the function to work something like:
  const double effStressInitialGuess = effStressT;
  double effStressTpdt =
    EffectiveStress::getEffStress(effStressInitialGuess,
				  effStressParams,
				  pylith::materials::PowerLaw3D::_effStressFunc,
				  pylith::materials::PowerLaw3D::_effStressFuncDFunc);

  // Compute Maxwell time
  properties[_PowerLaw3D::pidMaxwellTime] = 1.0e30;
  if (effStressTpdt != 0.0) {
    properties[_PowerLaw3D::pidMaxwellTime] =
      ((viscosityCoeff/effStressTpdt) ^ (powerLawExp - 1.0)) *
      (viscosityCoeff/mu);
    PetscLogFlops(5);
    } // if

    // Compute stresses from effective stress.
    const double effStressTau = (1.0 - alpha) * effStressT +
      alpha * effStressTpdt;
    const double gammaTau = 0.5 *
      ((effStressTau/viscosityCoeff)^(powerLawExp - 1.0))/visscosityCoeff;
    const double factor1 = 1.0/(ae + alpha * _dt * gammaTau);
    const double factor2 = timeFac * gammaTau;
    double devStressTpdt = 0.0;

    for (int iComp=0; iComp < tensorSize; ++iComp) {
      devStressTpdt = factor1 *
	(strainPPTpdt[iComp] - factor2 * devStressT[iComp] +
	 ae * devStressInitial[iComp]);
      stress[iComp] = devStressTpdt + diag[iComp] *
	(meanStressTpdt + meanStressInitial);
    } // for
    PetscLogFlops(14 + 8 * tensorSize);
  assert(0 != elasticConsts);
  assert(_PowerLaw3D::numElasticConsts == numElasticConsts);
  assert(0 != properties);
  assert(_totalPropsQuadPt == numProperties);
  assert(0 != totalStrain);
  assert(_PowerLaw3D::tensorSize == strainSize);
  assert(_PowerLaw3D::tensorSize == initialStateSize);
 
  const double mu = properties[_PowerLaw3D::pidMu];
  const double lambda = properties[_PowerLaw3D::pidLambda];
  const double maxwelltime = properties[_PowerLaw3D::pidMaxwellTime];

  const double mu2 = 2.0 * mu;
  const double bulkModulus = lambda + mu2/3.0;

  double dq = ViscoelasticMaxwell::computeVisStrain(_dt, maxwelltime);

  const double visFac = mu*dq/3.0;
  elasticConsts[ 0] = bulkModulus + 4.0*visFac; // C1111
  elasticConsts[ 1] = bulkModulus - 2.0*visFac; // C1122
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

  PetscLogFlops(10);
} // _calcElasticConstsViscoelastic

// ----------------------------------------------------------------------
// Update state variables.
void
pylith::materials::PowerLaw3D::_updatePropertiesElastic(
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
  assert(_PowerLaw3D::tensorSize == strainSize);

  const double maxwelltime = properties[_PowerLaw3D::pidMaxwellTime];

  const double e11 = totalStrain[0];
  const double e22 = totalStrain[1];
  const double e33 = totalStrain[2];

  const double traceStrainTpdt = e11 + e22 + e33;
  const double meanStrainTpdt = traceStrainTpdt/3.0;

  const double diag[] = { 1.0, 1.0, 1.0, 0.0, 0.0, 0.0 };

  for (int iComp=0; iComp < _PowerLaw3D::tensorSize; ++iComp) {
    properties[_PowerLaw3D::pidStrainT+iComp] = totalStrain[iComp];
    properties[_PowerLaw3D::pidVisStrainT+iComp] =
      totalStrain[iComp] - diag[iComp] * meanStrainTpdt;
  } // for
  PetscLogFlops(3 + 2 * _PowerLaw3D::tensorSize);

  _needNewJacobian = true;
} // _updatePropertiesElastic

// ----------------------------------------------------------------------
// Update state variables.
void
pylith::materials::PowerLaw3D::_updatePropertiesViscoelastic(
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
  assert(_PowerLaw3D::tensorSize == strainSize);
  assert(_PowerLaw3D::tensorSize == initialStateSize);

  const int tensorSize = _PowerLaw3D::tensorSize;

  pylith::materials::PowerLaw3D::_computeStateVars(properties,
							   numProperties,
							   totalStrain,
							   strainSize,
							   initialState,
							   initialStateSize);

  memcpy(&properties[_PowerLaw3D::pidVisStrainT],
	 &_visStrainT[0], 
	 tensorSize * sizeof(double));
  memcpy(&properties[_PowerLaw3D::pidStrainT],
	 &totalStrain[0], 
	 tensorSize * sizeof(double));

  _needNewJacobian = false;

} // _updatePropertiesViscoelastic


// End of file 
