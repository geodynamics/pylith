// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2015 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "MaxwellIsotropic3D.hh" // implementation of object methods

#include "ViscoelasticMaxwell.hh" // USES computeVisStrain
#include "Metadata.hh" // USES Metadata

#include "pylith/utils/array.hh" // USES scalar_array

#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include "petsc.h" // USES PetscLogFlops

#include <cassert> // USES assert()
#include <cstring> // USES memcpy()
#include <sstream> // USES std::ostringstream
#include <stdexcept> // USES std::runtime_error

// ----------------------------------------------------------------------
namespace pylith {
  namespace materials {
    namespace _MaxwellIsotropic3D{

      // Dimension of material.
      const int dimension = 3;

      /// Number of entries in stress/strain tensors.
      const int tensorSize = 6;

      /// Number of entries in derivative of elasticity matrix.
      const int numElasticConsts = 36;

      /// Number of physical properties.
      const int numProperties = 4;
      
      /// Physical properties.
      const Metadata::ParamDescription properties[] = {
	{ "density", 1, pylith::topology::FieldBase::SCALAR },
	{ "mu", 1, pylith::topology::FieldBase::SCALAR },
	{ "lambda", 1, pylith::topology::FieldBase::SCALAR },
	{ "maxwell_time", 1, pylith::topology::FieldBase::SCALAR },
      };
	
      // Values expected in properties spatial database
      const int numDBProperties = 4;
      const char* dbProperties[] = {"density", "vs", "vp" , "viscosity"};

      /// Number of state variables.
      const int numStateVars = 2;
      
      /// State variables.
      const Metadata::ParamDescription stateVars[] = {
	{ "total_strain", tensorSize, pylith::topology::FieldBase::TENSOR },
	{ "viscous_strain", tensorSize, pylith::topology::FieldBase::TENSOR },
      };

      // Values expected in state variables spatial database
      const int numDBStateVars = 2*tensorSize;
      const char* dbStateVars[] = {"total-strain-xx",
				   "total-strain-yy",
				   "total-strain-zz",
				   "total-strain-xy",
				   "total-strain-yz",
				   "total-strain-xz",
				   "viscous-strain-xx",
				   "viscous-strain-yy",
				   "viscous-strain-zz",
				   "viscous-strain-xy",
				   "viscous-strain-yz",
				   "viscous-strain-xz",
      };

    } // _MaxwellIsotropic3D
  } // materials
} // pylith

// Indices of physical properties
const int pylith::materials::MaxwellIsotropic3D::p_density = 0;

const int pylith::materials::MaxwellIsotropic3D::p_mu = 
  pylith::materials::MaxwellIsotropic3D::p_density + 1;

const int pylith::materials::MaxwellIsotropic3D::p_lambda = 
  pylith::materials::MaxwellIsotropic3D::p_mu + 1;

const int pylith::materials::MaxwellIsotropic3D::p_maxwellTime = 
  pylith::materials::MaxwellIsotropic3D::p_lambda + 1;

// Indices of database values (order must match dbProperties)
const int pylith::materials::MaxwellIsotropic3D::db_density = 0;

const int pylith::materials::MaxwellIsotropic3D::db_vs = 
  pylith::materials::MaxwellIsotropic3D::db_density + 1;

const int pylith::materials::MaxwellIsotropic3D::db_vp = 
  pylith::materials::MaxwellIsotropic3D::db_vs + 1;

const int pylith::materials::MaxwellIsotropic3D::db_viscosity = 
  pylith::materials::MaxwellIsotropic3D::db_vp + 1;

// Indices of state variables
const int pylith::materials::MaxwellIsotropic3D::s_totalStrain = 0;

const int pylith::materials::MaxwellIsotropic3D::s_viscousStrain = 
  pylith::materials::MaxwellIsotropic3D::s_totalStrain + 6;

// Indices of database values (order must match dbStateVars)
const int pylith::materials::MaxwellIsotropic3D::db_totalStrain = 0;

const int pylith::materials::MaxwellIsotropic3D::db_viscousStrain = 
  pylith::materials::MaxwellIsotropic3D::db_totalStrain + 6;

// ----------------------------------------------------------------------
// Default constructor.
pylith::materials::MaxwellIsotropic3D::MaxwellIsotropic3D(void) :
  ElasticMaterial(_MaxwellIsotropic3D::dimension,
		  _MaxwellIsotropic3D::tensorSize,
		  _MaxwellIsotropic3D::numElasticConsts,
		  Metadata(_MaxwellIsotropic3D::properties,
			   _MaxwellIsotropic3D::numProperties,
			   _MaxwellIsotropic3D::dbProperties,
			   _MaxwellIsotropic3D::numDBProperties,
			   _MaxwellIsotropic3D::stateVars,
			   _MaxwellIsotropic3D::numStateVars,
			   _MaxwellIsotropic3D::dbStateVars,
			   _MaxwellIsotropic3D::numDBStateVars)),
  _calcElasticConstsFn(0),
  _calcStressFn(0),
  _updateStateVarsFn(0)
{ // constructor
  useElasticBehavior(false);
  _viscousStrain.resize(_tensorSize);
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::materials::MaxwellIsotropic3D::~MaxwellIsotropic3D(void)
{ // destructor
} // destructor

// ----------------------------------------------------------------------
// Set whether elastic or inelastic constitutive relations are used.
void
pylith::materials::MaxwellIsotropic3D::useElasticBehavior(const bool flag)
{ // useElasticBehavior
  if (flag) {
    _calcStressFn = 
      &pylith::materials::MaxwellIsotropic3D::_calcStressElastic;
    _calcElasticConstsFn = 
      &pylith::materials::MaxwellIsotropic3D::_calcElasticConstsElastic;
    _updateStateVarsFn = 
      &pylith::materials::MaxwellIsotropic3D::_updateStateVarsElastic;

  } else {
    _calcStressFn = 
      &pylith::materials::MaxwellIsotropic3D::_calcStressViscoelastic;
    _calcElasticConstsFn = 
      &pylith::materials::MaxwellIsotropic3D::_calcElasticConstsViscoelastic;
    _updateStateVarsFn = 
      &pylith::materials::MaxwellIsotropic3D::_updateStateVarsViscoelastic;
  } // if/else
} // useElasticBehavior

// ----------------------------------------------------------------------
// Compute properties from values in spatial database.
void
pylith::materials::MaxwellIsotropic3D::_dbToProperties(
					    PylithScalar* const propValues,
					    const scalar_array& dbValues)
{ // _dbToProperties
  assert(0 != propValues);
  const int numDBValues = dbValues.size();
  assert(_MaxwellIsotropic3D::numDBProperties == numDBValues);

  const PylithScalar density = dbValues[db_density];
  const PylithScalar vs = dbValues[db_vs];
  const PylithScalar vp = dbValues[db_vp];
  const PylithScalar viscosity = dbValues[db_viscosity];
 
  if (density <= 0.0 || vs <= 0.0 || vp <= 0.0 || viscosity <= 0.0) {
    std::ostringstream msg;
    msg << "Spatial database returned nonpositive value for physical "
	<< "properties.\n"
	<< "density: " << density << "\n"
	<< "vp: " << vp << "\n"
	<< "vs: " << vs << "\n"
	<< "viscosity: " << viscosity << "\n";
    throw std::runtime_error(msg.str());
  } // if

  const PylithScalar mu = density * vs*vs;
  const PylithScalar lambda = density * vp*vp - 2.0*mu;

  if (lambda <= 0.0) {
    std::ostringstream msg;
    msg << "Attempted to set Lame's constant lambda to nonpositive value.\n"
	<< "density: " << density << "\n"
	<< "vp: " << vp << "\n"
	<< "vs: " << vs << "\n";
    throw std::runtime_error(msg.str());
  } // if
  assert(mu > 0);

  const PylithScalar maxwellTime = viscosity / mu;
  assert(maxwellTime > 0.0);

  propValues[p_density] = density;
  propValues[p_mu] = mu;
  propValues[p_lambda] = lambda;
  propValues[p_maxwellTime] = maxwellTime;

  PetscLogFlops(7);
} // _dbToProperties

// ----------------------------------------------------------------------
// Nondimensionalize properties.
void
pylith::materials::MaxwellIsotropic3D::_nondimProperties(PylithScalar* const values,
							 const int nvalues) const
{ // _nondimProperties
  assert(0 != _normalizer);
  assert(0 != values);
  assert(nvalues == _numPropsQuadPt);

  const PylithScalar densityScale = _normalizer->densityScale();
  const PylithScalar pressureScale = _normalizer->pressureScale();
  const PylithScalar timeScale = _normalizer->timeScale();
  values[p_density] = 
    _normalizer->nondimensionalize(values[p_density], densityScale);
  values[p_mu] = 
    _normalizer->nondimensionalize(values[p_mu], pressureScale);
  values[p_lambda] = 
    _normalizer->nondimensionalize(values[p_lambda], pressureScale);
  values[p_maxwellTime] = 
    _normalizer->nondimensionalize(values[p_maxwellTime], timeScale);

  PetscLogFlops(4);
} // _nondimProperties

// ----------------------------------------------------------------------
// Dimensionalize properties.
void
pylith::materials::MaxwellIsotropic3D::_dimProperties(PylithScalar* const values,
						      const int nvalues) const
{ // _dimProperties
  assert(0 != _normalizer);
  assert(0 != values);
  assert(nvalues == _numPropsQuadPt);

  const PylithScalar densityScale = _normalizer->densityScale();
  const PylithScalar pressureScale = _normalizer->pressureScale();
  const PylithScalar timeScale = _normalizer->timeScale();
  values[p_density] = 
    _normalizer->dimensionalize(values[p_density], densityScale);
  values[p_mu] = 
    _normalizer->dimensionalize(values[p_mu], pressureScale);
  values[p_lambda] = 
    _normalizer->dimensionalize(values[p_lambda], pressureScale);
  values[p_maxwellTime] = 
    _normalizer->dimensionalize(values[p_maxwellTime], timeScale);

  PetscLogFlops(4);
} // _dimProperties

// ----------------------------------------------------------------------
// Compute initial state variables from values in spatial database.
void
pylith::materials::MaxwellIsotropic3D::_dbToStateVars(
					PylithScalar* const stateValues,
					const scalar_array& dbValues)
{ // _dbToStateVars
  assert(0 != stateValues);
  const int numDBValues = dbValues.size();
  assert(_MaxwellIsotropic3D::numDBStateVars == numDBValues);

  const int totalSize = 2 * _tensorSize;
  assert(totalSize == _numVarsQuadPt);
  assert(totalSize == numDBValues);
  memcpy(stateValues, &dbValues[0], totalSize*sizeof(PylithScalar));

  PetscLogFlops(0);
} // _dbToStateVars

// ----------------------------------------------------------------------
// Compute density at location from properties.
void
pylith::materials::MaxwellIsotropic3D::_calcDensity(PylithScalar* const density,
						    const PylithScalar* properties,
						    const int numProperties,
						    const PylithScalar* stateVars,
						    const int numStateVars)
{ // _calcDensity
  assert(0 != density);
  assert(0 != properties);
  assert(_numPropsQuadPt == numProperties);

  density[0] = properties[p_density];
} // _calcDensity

// ----------------------------------------------------------------------
// Compute stress tensor at location from properties as an elastic
// material.
void
pylith::materials::MaxwellIsotropic3D::_calcStressElastic(
					     PylithScalar* const stress,
					     const int stressSize,
					     const PylithScalar* properties,
					     const int numProperties,
					     const PylithScalar* stateVars,
					     const int numStateVars,
					     const PylithScalar* totalStrain,
					     const int strainSize,
					     const PylithScalar* initialStress,
					     const int initialStressSize,
					     const PylithScalar* initialStrain,
					     const int initialStrainSize,
					     const bool computeStateVars)
{ // _calcStressElastic
  assert(0 != stress);
  assert(_MaxwellIsotropic3D::tensorSize == stressSize);
  assert(0 != properties);
  assert(_numPropsQuadPt == numProperties);
  assert(0 != stateVars);
  assert(_numVarsQuadPt == numStateVars);
  assert(0 != totalStrain);
  assert(_MaxwellIsotropic3D::tensorSize == strainSize);
  assert(0 != initialStress);
  assert(_MaxwellIsotropic3D::tensorSize == initialStressSize);
  assert(0 != initialStrain);
  assert(_MaxwellIsotropic3D::tensorSize == initialStrainSize);

  const PylithScalar mu = properties[p_mu];
  const PylithScalar lambda = properties[p_lambda];
  const PylithScalar mu2 = 2.0 * mu;

  // :TODO: Need to consider initial state variables????
  const PylithScalar e11 = totalStrain[0] - initialStrain[0];
  const PylithScalar e22 = totalStrain[1] - initialStrain[1];
  const PylithScalar e33 = totalStrain[2] - initialStrain[2];
  const PylithScalar e12 = totalStrain[3] - initialStrain[3];
  const PylithScalar e23 = totalStrain[4] - initialStrain[4];
  const PylithScalar e13 = totalStrain[5] - initialStrain[5];
  
  const PylithScalar s123 = lambda * (e11 + e22 + e33);

  stress[0] = s123 + mu2 * e11 + initialStress[0];
  stress[1] = s123 + mu2 * e22 + initialStress[1];
  stress[2] = s123 + mu2 * e33 + initialStress[2];
  stress[3] = mu2 * e12 + initialStress[3];
  stress[4] = mu2 * e23 + initialStress[4];
  stress[5] = mu2 * e13 + initialStress[5];

  PetscLogFlops(25);
} // _calcStressElastic

// ----------------------------------------------------------------------
// Compute stress tensor at location from properties as a viscoelastic
// material.
void
pylith::materials::MaxwellIsotropic3D::_calcStressViscoelastic(
					     PylithScalar* const stress,
					     const int stressSize,
					     const PylithScalar* properties,
					     const int numProperties,
					     const PylithScalar* stateVars,
					     const int numStateVars,
					     const PylithScalar* totalStrain,
					     const int strainSize,
					     const PylithScalar* initialStress,
					     const int initialStressSize,
					     const PylithScalar* initialStrain,
					     const int initialStrainSize,
					     const bool computeStateVars)
{ // _calcStressViscoelastic
  assert(0 != stress);
  assert(_MaxwellIsotropic3D::tensorSize == stressSize);
  assert(0 != properties);
  assert(_numPropsQuadPt == numProperties);
  assert(0 != stateVars);
  assert(_numVarsQuadPt == numStateVars);
  assert(0 != totalStrain);
  assert(_MaxwellIsotropic3D::tensorSize == strainSize);
  assert(0 != initialStress);
  assert(_MaxwellIsotropic3D::tensorSize == initialStressSize);
  assert(0 != initialStrain);
  assert(_MaxwellIsotropic3D::tensorSize == initialStrainSize);

  const int tensorSize = _MaxwellIsotropic3D::tensorSize;

  const PylithScalar mu = properties[p_mu];
  const PylithScalar lambda = properties[p_lambda];

  const PylithScalar mu2 = 2.0 * mu;
  const PylithScalar bulkModulus = lambda + mu2 / 3.0;

  // Initial stress and strain values
  const PylithScalar meanStrainInitial = (initialStrain[0] +
				    initialStrain[1] +
				    initialStrain[2]) / 3.0;
  const PylithScalar meanStressInitial = (initialStress[0] +
				    initialStress[1] +
				    initialStress[2]) / 3.0;
  const PylithScalar devStrainInitial[] = {initialStrain[0] - meanStrainInitial,
				     initialStrain[1] - meanStrainInitial,
				     initialStrain[2] - meanStrainInitial,
				     initialStrain[3],
				     initialStrain[4],
				     initialStrain[5]};
  const PylithScalar devStressInitial[] = {initialStress[0] - meanStressInitial,
				     initialStress[1] - meanStressInitial,
				     initialStress[2] - meanStressInitial,
				     initialStress[3],
				     initialStress[4],
				     initialStress[5]};

  // :TODO: Need to determine how to incorporate state variables
  const PylithScalar meanStrainTpdt = (totalStrain[0] +
				 totalStrain[1] +
				 totalStrain[2]) / 3.0;
  
  const PylithScalar meanStressTpdt = 3.0 * bulkModulus *
    (meanStrainTpdt - meanStrainInitial) + meanStressInitial;

  const PylithScalar diag[] = { 1.0, 1.0, 1.0, 0.0, 0.0, 0.0 };

  // Get viscous strains
  if (computeStateVars)
    _computeStateVars(stateVars, numStateVars,
		      properties, numProperties,
		      totalStrain, strainSize,
		      initialStress, initialStressSize,
		      initialStrain, initialStrainSize);
  else
    memcpy(&_viscousStrain[0], &stateVars[s_viscousStrain],
	   tensorSize*sizeof(PylithScalar));

  // Compute new stresses
  PylithScalar devStressTpdt = 0.0;

  for (int iComp=0; iComp < tensorSize; ++iComp) {
    devStressTpdt = mu2 * (_viscousStrain[iComp] - devStrainInitial[iComp]) +
      devStressInitial[iComp];

    stress[iComp] = diag[iComp] * meanStressTpdt + devStressTpdt;
  } // for

  PetscLogFlops(22 + 5 * tensorSize);
} // _calcStressViscoelastic

// ----------------------------------------------------------------------
// Compute derivative of elasticity matrix at location from properties.
void
pylith::materials::MaxwellIsotropic3D::_calcElasticConstsElastic(
					PylithScalar* const elasticConsts,
					const int numElasticConsts,
					const PylithScalar* properties,
					const int numProperties,
					const PylithScalar* stateVars,
					const int numStateVars,
					const PylithScalar* totalStrain,
					const int strainSize,
					const PylithScalar* initialStress,
					const int initialStressSize,
					const PylithScalar* initialStrain,
					const int initialStrainSize)
{ // _calcElasticConstsElastic
  assert(0 != elasticConsts);
  assert(_MaxwellIsotropic3D::numElasticConsts == numElasticConsts);
  assert(0 != properties);
  assert(_numPropsQuadPt == numProperties);
  assert(0 != stateVars);
  assert(_numVarsQuadPt == numStateVars);
  assert(0 != totalStrain);
  assert(_MaxwellIsotropic3D::tensorSize == strainSize);
  assert(0 != initialStress);
  assert(_MaxwellIsotropic3D::tensorSize == initialStressSize);
  assert(0 != initialStrain);
  assert(_MaxwellIsotropic3D::tensorSize == initialStrainSize);
 
  const PylithScalar mu = properties[p_mu];
  const PylithScalar lambda = properties[p_lambda];

  const PylithScalar mu2 = 2.0 * mu;
  const PylithScalar lambda2mu = lambda + mu2;

  elasticConsts[ 0] = lambda2mu; // C1111
  elasticConsts[ 1] = lambda; // C1122
  elasticConsts[ 2] = lambda; // C1133
  elasticConsts[ 3] = 0; // C1112
  elasticConsts[ 4] = 0; // C1123
  elasticConsts[ 5] = 0; // C1113
  elasticConsts[ 6] = lambda; // C2211
  elasticConsts[ 7] = lambda2mu; // C2222
  elasticConsts[ 8] = lambda; // C2233
  elasticConsts[ 9] = 0; // C2212
  elasticConsts[10] = 0; // C2223
  elasticConsts[11] = 0; // C2213
  elasticConsts[12] = lambda; // C3311
  elasticConsts[13] = lambda; // C3322
  elasticConsts[14] = lambda2mu; // C3333
  elasticConsts[15] = 0; // C3312
  elasticConsts[16] = 0; // C3323
  elasticConsts[17] = 0; // C3313
  elasticConsts[18] = 0; // C1211
  elasticConsts[19] = 0; // C1222
  elasticConsts[20] = 0; // C1233
  elasticConsts[21] = mu2; // C1212
  elasticConsts[22] = 0; // C1223
  elasticConsts[23] = 0; // C1213
  elasticConsts[24] = 0; // C2311
  elasticConsts[25] = 0; // C2322
  elasticConsts[26] = 0; // C2333
  elasticConsts[27] = 0; // C2312
  elasticConsts[28] = mu2; // C2323
  elasticConsts[29] = 0; // C2313
  elasticConsts[30] = 0; // C1311
  elasticConsts[31] = 0; // C1322
  elasticConsts[32] = 0; // C1333
  elasticConsts[33] = 0; // C1312
  elasticConsts[34] = 0; // C1323
  elasticConsts[35] = mu2; // C1313

  PetscLogFlops(2);
} // _calcElasticConstsElastic

// ----------------------------------------------------------------------
// Compute derivative of elasticity matrix at location from properties
// as an elastic material.
void
pylith::materials::MaxwellIsotropic3D::_calcElasticConstsViscoelastic(
					PylithScalar* const elasticConsts,
					const int numElasticConsts,
					const PylithScalar* properties,
					const int numProperties,
					const PylithScalar* stateVars,
					const int numStateVars,
					const PylithScalar* totalStrain,
					const int strainSize,
					const PylithScalar* initialStress,
					const int initialStressSize,
					const PylithScalar* initialStrain,
					const int initialStrainSize)
{ // _calcElasticConstsViscoelastic
  assert(0 != elasticConsts);
  assert(_MaxwellIsotropic3D::numElasticConsts == numElasticConsts);
  assert(0 != properties);
  assert(_numPropsQuadPt == numProperties);
  assert(0 != stateVars);
  assert(_numVarsQuadPt == numStateVars);
  assert(0 != totalStrain);
  assert(_MaxwellIsotropic3D::tensorSize == strainSize);
  assert(0 != initialStress);
  assert(_MaxwellIsotropic3D::tensorSize == initialStressSize);
  assert(0 != initialStrain);
  assert(_MaxwellIsotropic3D::tensorSize == initialStrainSize);

  const PylithScalar mu = properties[p_mu];
  const PylithScalar lambda = properties[p_lambda];
  const PylithScalar maxwellTime = properties[p_maxwellTime];

  const PylithScalar mu2 = 2.0 * mu;
  const PylithScalar bulkModulus = lambda + mu2 / 3.0;

  PylithScalar dq = ViscoelasticMaxwell::viscousStrainParam(_dt, maxwellTime);

  const PylithScalar visFac = mu * dq / 3.0;

  elasticConsts[ 0] = bulkModulus + 4.0 * visFac; // C1111
  elasticConsts[ 1] = bulkModulus - 2.0 * visFac; // C1122
  elasticConsts[ 2] = elasticConsts[1]; // C1133
  elasticConsts[ 3] = 0; // C1112
  elasticConsts[ 4] = 0; // C1123
  elasticConsts[ 5] = 0; // C1113
  elasticConsts[ 6] = elasticConsts[1]; // C2211
  elasticConsts[ 7] = elasticConsts[0]; // C2222
  elasticConsts[ 8] = elasticConsts[1]; // C2233
  elasticConsts[ 9] = 0; // C2212
  elasticConsts[10] = 0; // C2223
  elasticConsts[11] = 0; // C2213
  elasticConsts[12] = elasticConsts[1]; // C3311
  elasticConsts[13] = elasticConsts[1]; // C3322
  elasticConsts[14] = elasticConsts[0]; // C3333
  elasticConsts[15] = 0; // C3312
  elasticConsts[16] = 0; // C3323
  elasticConsts[17] = 0; // C3313
  elasticConsts[18] = 0; // C1211
  elasticConsts[19] = 0; // C1222
  elasticConsts[20] = 0; // C1233
  elasticConsts[21] = 6.0 * visFac; // C1212
  elasticConsts[22] = 0; // C1223
  elasticConsts[23] = 0; // C1213
  elasticConsts[24] = 0; // C2311
  elasticConsts[25] = 0; // C2322
  elasticConsts[26] = 0; // C2333
  elasticConsts[27] = 0; // C2312
  elasticConsts[28] = elasticConsts[21]; // C2323
  elasticConsts[29] = 0; // C2313
  elasticConsts[30] = 0; // C1311
  elasticConsts[31] = 0; // C1322
  elasticConsts[32] = 0; // C1333
  elasticConsts[33] = 0; // C1312
  elasticConsts[34] = 0; // C1323
  elasticConsts[35] = elasticConsts[21]; // C1313

  PetscLogFlops(10);
} // _calcElasticConstsViscoelastic

// ----------------------------------------------------------------------
// Update state variables as an elastic material.
void
pylith::materials::MaxwellIsotropic3D::_updateStateVarsElastic(
					    PylithScalar* const stateVars,
					    const int numStateVars,
					    const PylithScalar* properties,
					    const int numProperties,
					    const PylithScalar* totalStrain,
					    const int strainSize,
					    const PylithScalar* initialStress,
					    const int initialStressSize,
					    const PylithScalar* initialStrain,
					    const int initialStrainSize)
{ // _updateStateVarsElastic
  assert(0 != stateVars);
  assert(_numVarsQuadPt == numStateVars);
  assert(0 != properties);
  assert(_numPropsQuadPt == numProperties);
  assert(0 != totalStrain);
  assert(_MaxwellIsotropic3D::tensorSize == strainSize);
  assert(0 != initialStress);
  assert(_MaxwellIsotropic3D::tensorSize == initialStressSize);
  assert(0 != initialStrain);
  assert(_MaxwellIsotropic3D::tensorSize == initialStrainSize);

  const int tensorSize = _tensorSize;

  const PylithScalar strainTpdt[] = {totalStrain[0] - initialStrain[0],
			       totalStrain[1] - initialStrain[1],
			       totalStrain[2] - initialStrain[2],
			       totalStrain[3] - initialStrain[3],
			       totalStrain[4] - initialStrain[4],
			       totalStrain[5] - initialStrain[5]};

  const PylithScalar meanStrainTpdt =
    (strainTpdt[0] + strainTpdt[1] + strainTpdt[2]) / 3.0;

  const PylithScalar diag[] = { 1.0, 1.0, 1.0, 0.0, 0.0, 0.0 };

  // Initialize all viscous strains to deviatoric elastic strains.
  for (int iComp=0; iComp < tensorSize; ++iComp) {
    stateVars[s_totalStrain+iComp] = totalStrain[iComp];
    stateVars[s_viscousStrain+iComp] =
      strainTpdt[iComp] - diag[iComp] * meanStrainTpdt;
  } // for
  PetscLogFlops(9 + 2 * _tensorSize);

  _needNewJacobian = true;
} // _updateStateVarsElastic

// ----------------------------------------------------------------------
// Update state variables as a viscoelastic material.
void
pylith::materials::MaxwellIsotropic3D::_updateStateVarsViscoelastic(
					    PylithScalar* const stateVars,
					    const int numStateVars,
					    const PylithScalar* properties,
					    const int numProperties,
					    const PylithScalar* totalStrain,
					    const int strainSize,
					    const PylithScalar* initialStress,
					    const int initialStressSize,
					    const PylithScalar* initialStrain,
					    const int initialStrainSize)
{ // _updateStateVarsViscoelastic
  assert(0 != stateVars);
  assert(_numVarsQuadPt == numStateVars);
  assert(0 != properties);
  assert(_numPropsQuadPt == numProperties);
  assert(0 != totalStrain);
  assert(_MaxwellIsotropic3D::tensorSize == strainSize);
  assert(0 != initialStress);
  assert(_MaxwellIsotropic3D::tensorSize == initialStressSize);
  assert(0 != initialStrain);
  assert(_MaxwellIsotropic3D::tensorSize == initialStrainSize);

  const int tensorSize = _tensorSize;

  _computeStateVars(stateVars, numStateVars,
		    properties, numProperties,
		    totalStrain, strainSize,
		    initialStress, initialStressSize,
		    initialStrain, initialStrainSize);

  memcpy(&stateVars[s_totalStrain], totalStrain, tensorSize*sizeof(PylithScalar));

  memcpy(&stateVars[s_viscousStrain], &_viscousStrain[0], 
	 tensorSize*sizeof(PylithScalar));

  _needNewJacobian = false;
} // _updateStateVarsViscoelastic

// ----------------------------------------------------------------------
// Get stable time step for implicit time integration.
PylithScalar
pylith::materials::MaxwellIsotropic3D::_stableTimeStepImplicit(
					   const PylithScalar* properties,
					   const int numProperties,
					   const PylithScalar* stateVars,
					   const int numStateVars) const
{ // _stableTimeStepImplicit
  assert(0 != properties);
  assert(_numPropsQuadPt == numProperties);
  assert(0 != stateVars);
  assert(_numVarsQuadPt == numStateVars);

  const PylithScalar maxwellTime = properties[p_maxwellTime];
  const PylithScalar dtStable = 0.2 * maxwellTime;

  return dtStable;
} // _stableTimeStepImplicit

// ----------------------------------------------------------------------
// Get stable time step for explicit time integration.
PylithScalar
pylith::materials::MaxwellIsotropic3D::_stableTimeStepExplicit(const PylithScalar* properties,
							       const int numProperties,
							       const PylithScalar* stateVars,
							       const int numStateVars,
							       const double minCellWidth) const
{ // _stableTimeStepExplicit
  assert(properties);
  assert(_numPropsQuadPt == numProperties);
 
  const PylithScalar mu = properties[p_mu];
  const PylithScalar lambda = properties[p_lambda];
  const PylithScalar density = properties[p_density];

  assert(density > 0.0);
  const PylithScalar vp = sqrt((lambda + 2*mu) / density);

  const PylithScalar dtStable = minCellWidth / vp;
  return dtStable;
} // _stableTimeStepExplicit


// ----------------------------------------------------------------------
// Compute viscous strain for current time step.
void
pylith::materials::MaxwellIsotropic3D::_computeStateVars(
					       const PylithScalar* stateVars,
					       const int numStateVars,
					       const PylithScalar* properties,
					       const int numProperties,
					       const PylithScalar* totalStrain,
					       const int strainSize,
					       const PylithScalar* initialStress,
					       const int initialStressSize,
					       const PylithScalar* initialStrain,
					       const int initialStrainSize)
{ // _computeStateVars
  assert(0 != stateVars);
  assert(_numVarsQuadPt == numStateVars);
  assert(0 != properties);
  assert(_numPropsQuadPt == numProperties);
  assert(0 != totalStrain);
  assert(_MaxwellIsotropic3D::tensorSize == strainSize);
  assert(0 != initialStress);
  assert(_MaxwellIsotropic3D::tensorSize == initialStressSize);
  assert(0 != initialStrain);
  assert(_MaxwellIsotropic3D::tensorSize == initialStrainSize);

  const int tensorSize = _tensorSize;
  const PylithScalar maxwellTime = properties[p_maxwellTime];

  // :TODO: Need to account for initial values for state variables
  // and the initial strain??

  const PylithScalar e11 = totalStrain[0];
  const PylithScalar e22 = totalStrain[1];
  const PylithScalar e33 = totalStrain[2];
  
  const PylithScalar meanStrainTpdt = (e11 + e22 + e33) / 3.0;

  const PylithScalar diag[] = { 1.0, 1.0, 1.0, 0.0, 0.0, 0.0 };

  const PylithScalar meanStrainT =
    ( stateVars[s_totalStrain+0] +
      stateVars[s_totalStrain+1] +
      stateVars[s_totalStrain+2] ) / 3.0;
  
  // Time integration.
  PylithScalar dq = ViscoelasticMaxwell::viscousStrainParam(_dt, maxwellTime);
  const PylithScalar expFac = exp(-_dt/maxwellTime);

  PylithScalar devStrainTpdt = 0.0;
  PylithScalar devStrainT = 0.0;

  for (int iComp=0; iComp < tensorSize; ++iComp) {
    devStrainTpdt = totalStrain[iComp] - diag[iComp] * meanStrainTpdt;
    devStrainT = stateVars[s_totalStrain+iComp] - diag[iComp] * meanStrainT;
    _viscousStrain[iComp] = expFac * stateVars[s_viscousStrain+iComp] +
      dq * (devStrainTpdt - devStrainT);
  } // for

  PetscLogFlops(9 + 7 * tensorSize);
} // _computeStateVars


// End of file 
