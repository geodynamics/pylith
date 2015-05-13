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

#include "GenMaxwellIsotropic3D.hh" // implementation of object methods

#include "ViscoelasticMaxwell.hh" // USES computeVisStrain
#include "Metadata.hh" // USES Metadata

#include "pylith/utils/array.hh" // USES scalar_array
#include "pylith/utils/constdefs.h" // USES MAXSCALAR

#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include "petsc.h" // USES PetscLogFlops

#include <cassert> // USES assert()
#include <sstream> // USES std::ostringstream
#include <stdexcept> // USES std::runtime_error

// ----------------------------------------------------------------------
namespace pylith {
  namespace materials {
    namespace _GenMaxwellIsotropic3D{

      /// Number of Maxwell models in parallel.
      const int numMaxwellModels = 3;

      // Dimension of material.
      const int dimension = 3;

      /// Number of entries in stress/strain tensors.
      const int tensorSize = 6;

      /// Number of entries in derivative of elasticity matrix.
      const int numElasticConsts = 36;

      /// Number of physical properties.
      const int numProperties = 5;
      
      /// Physical properties.
      const Metadata::ParamDescription properties[numProperties] = {
	{ "density", 1, pylith::topology::FieldBase::SCALAR },
	{ "mu", 1, pylith::topology::FieldBase::SCALAR },
	{ "lambda", 1, pylith::topology::FieldBase::SCALAR },
	{ "shear_ratio", numMaxwellModels, pylith::topology::FieldBase::OTHER },
	{ "maxwell_time", numMaxwellModels, pylith::topology::FieldBase::OTHER },
      };
      // Values expected in properties spatial database.  :KLUDGE: Not
      // generalized over number of models.
      const int numDBProperties = 3 + 2*numMaxwellModels;
      const char* dbProperties[numDBProperties] = {
	"density", "vs", "vp",
	"shear-ratio-1",
	"shear-ratio-2",
	"shear-ratio-3",
	"viscosity-1",
	"viscosity-2",
	"viscosity-3",
      };
      
      /// Number of state variables.
      const int numStateVars = 1 + numMaxwellModels;
      
      /// State variables. :KLUDGE: Not generalized over number of models.
      const Metadata::ParamDescription stateVars[numStateVars] = {
	{ "total_strain", tensorSize, pylith::topology::FieldBase::TENSOR },
	{ "viscous_strain_1", tensorSize, pylith::topology::FieldBase::TENSOR },
	{ "viscous_strain_2", tensorSize, pylith::topology::FieldBase::TENSOR },
	{ "viscous_strain_3", tensorSize, pylith::topology::FieldBase::TENSOR },
      };

      // Values expected in state variables spatial database
      const int numDBStateVars = tensorSize + numMaxwellModels * tensorSize;
      const char* dbStateVars[numDBStateVars] = {"total-strain-xx",
						 "total-strain-yy",
						 "total-strain-zz",
						 "total-strain-xy",
						 "total-strain-yz",
						 "total-strain-xz",
						 "viscous-strain-1-xx",
						 "viscous-strain-1-yy",
						 "viscous-strain-1-zz",
						 "viscous-strain-1-xy",
						 "viscous-strain-1-yz",
						 "viscous-strain-1-xz",
						 "viscous-strain-2-xx",
						 "viscous-strain-2-yy",
						 "viscous-strain-2-zz",
						 "viscous-strain-2-xy",
						 "viscous-strain-2-yz",
						 "viscous-strain-2-xz",
						 "viscous-strain-3-xx",
						 "viscous-strain-3-yy",
						 "viscous-strain-3-zz",
						 "viscous-strain-3-xy",
						 "viscous-strain-3-yz",
						 "viscous-strain-3-xz",
      };

    } // _GenMaxwellIsotropic3D
  } // materials
} // pylith

// Indices of physical properties
const int pylith::materials::GenMaxwellIsotropic3D::p_density = 0;

const int pylith::materials::GenMaxwellIsotropic3D::p_muEff =
  pylith::materials::GenMaxwellIsotropic3D::p_density + 1;

const int pylith::materials::GenMaxwellIsotropic3D::p_lambdaEff =
  pylith::materials::GenMaxwellIsotropic3D::p_muEff + 1;

const int pylith::materials::GenMaxwellIsotropic3D::p_shearRatio =
  pylith::materials::GenMaxwellIsotropic3D::p_lambdaEff + 1;

const int pylith::materials::GenMaxwellIsotropic3D::p_maxwellTime =
  pylith::materials::GenMaxwellIsotropic3D::p_shearRatio +
  pylith::materials::_GenMaxwellIsotropic3D::numMaxwellModels;

// Indices of database values (order must match dbProperties)
const int pylith::materials::GenMaxwellIsotropic3D::db_density = 0;

const int pylith::materials::GenMaxwellIsotropic3D::db_vs =
  pylith::materials::GenMaxwellIsotropic3D::db_density + 1;

const int pylith::materials::GenMaxwellIsotropic3D::db_vp =
  pylith::materials::GenMaxwellIsotropic3D::db_vs + 1;

const int pylith::materials::GenMaxwellIsotropic3D::db_shearRatio =
  pylith::materials::GenMaxwellIsotropic3D::db_vp + 1;

const int pylith::materials::GenMaxwellIsotropic3D::db_viscosity =
  pylith::materials::GenMaxwellIsotropic3D::db_shearRatio + 
  pylith::materials::_GenMaxwellIsotropic3D::numMaxwellModels;

// Indices of state variables
const int pylith::materials::GenMaxwellIsotropic3D::s_totalStrain = 0;

const int pylith::materials::GenMaxwellIsotropic3D::s_viscousStrain1 = 
  pylith::materials::GenMaxwellIsotropic3D::s_totalStrain +
  pylith::materials::_GenMaxwellIsotropic3D::tensorSize;

const int pylith::materials::GenMaxwellIsotropic3D::s_viscousStrain2 = 
  pylith::materials::GenMaxwellIsotropic3D::s_viscousStrain1 +
  pylith::materials::_GenMaxwellIsotropic3D::tensorSize;

const int pylith::materials::GenMaxwellIsotropic3D::s_viscousStrain3 = 
  pylith::materials::GenMaxwellIsotropic3D::s_viscousStrain2 +
  pylith::materials::_GenMaxwellIsotropic3D::tensorSize;

// Indices of database values (order must match dbStateVars)
const int pylith::materials::GenMaxwellIsotropic3D::db_totalStrain = 0;

const int pylith::materials::GenMaxwellIsotropic3D::db_viscousStrain1 =
  pylith::materials::GenMaxwellIsotropic3D::db_totalStrain +
  pylith::materials::_GenMaxwellIsotropic3D::tensorSize;

const int pylith::materials::GenMaxwellIsotropic3D::db_viscousStrain2 =
  pylith::materials::GenMaxwellIsotropic3D::db_viscousStrain1 +
  pylith::materials::_GenMaxwellIsotropic3D::tensorSize;

const int pylith::materials::GenMaxwellIsotropic3D::db_viscousStrain3 =
  pylith::materials::GenMaxwellIsotropic3D::db_viscousStrain2 +
  pylith::materials::_GenMaxwellIsotropic3D::tensorSize;

// ----------------------------------------------------------------------
// Default constructor.
pylith::materials::GenMaxwellIsotropic3D::GenMaxwellIsotropic3D(void) :
  ElasticMaterial(_GenMaxwellIsotropic3D::dimension,
		  _GenMaxwellIsotropic3D::tensorSize,
		  _GenMaxwellIsotropic3D::numElasticConsts,
		  Metadata(_GenMaxwellIsotropic3D::properties,
			   _GenMaxwellIsotropic3D::numProperties,
			   _GenMaxwellIsotropic3D::dbProperties,
			   _GenMaxwellIsotropic3D::numDBProperties,
			   _GenMaxwellIsotropic3D::stateVars,
			   _GenMaxwellIsotropic3D::numStateVars,
			   _GenMaxwellIsotropic3D::dbStateVars,
			   _GenMaxwellIsotropic3D::numDBStateVars)),
  _calcElasticConstsFn(0),
  _calcStressFn(0),
  _updateStateVarsFn(0)  
{ // constructor
  useElasticBehavior(false);
  _viscousStrain.resize(_GenMaxwellIsotropic3D::numMaxwellModels*_tensorSize);
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::materials::GenMaxwellIsotropic3D::~GenMaxwellIsotropic3D(void)
{ // destructor
} // destructor

// ----------------------------------------------------------------------
// Set whether elastic or inelastic constitutive relations are used.
void
pylith::materials::GenMaxwellIsotropic3D::useElasticBehavior(const bool flag)
{ // useElasticBehavior
  if (flag) {
    _calcStressFn = 
      &pylith::materials::GenMaxwellIsotropic3D::_calcStressElastic;
    _calcElasticConstsFn = 
      &pylith::materials::GenMaxwellIsotropic3D::_calcElasticConstsElastic;
    _updateStateVarsFn = 
      &pylith::materials::GenMaxwellIsotropic3D::_updateStateVarsElastic;

  } else {
    _calcStressFn = 
      &pylith::materials::GenMaxwellIsotropic3D::_calcStressViscoelastic;
    _calcElasticConstsFn = 
      &pylith::materials::GenMaxwellIsotropic3D::_calcElasticConstsViscoelastic;
    _updateStateVarsFn = 
      &pylith::materials::GenMaxwellIsotropic3D::_updateStateVarsViscoelastic;
  } // if/else
} // useElasticBehavior

// ----------------------------------------------------------------------
// Compute parameters from values in spatial database.
void
pylith::materials::GenMaxwellIsotropic3D::_dbToProperties(
					    PylithScalar* const propValues,
					    const scalar_array& dbValues)
{ // _dbToProperties
  assert(0 != propValues);
  const int numDBValues = dbValues.size();
  assert(_GenMaxwellIsotropic3D::numDBProperties == numDBValues);

  const PylithScalar density = dbValues[db_density];
  const PylithScalar vs = dbValues[db_vs];
  const PylithScalar vp = dbValues[db_vp];
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

  propValues[p_density] = density;
  propValues[p_muEff] = mu;
  propValues[p_lambdaEff] = lambda;

  PylithScalar visFrac = 0.0;
  for (int imodel = 0; imodel < numMaxwellModels; ++imodel) 
    visFrac += dbValues[db_shearRatio + imodel];
  if (visFrac > 1.0) {
    std::ostringstream msg;
    msg << "Shear modulus ratios sum to a value greater than 1.0 for\n"
	<< "Generalized Maxwell model.\n"
	<< "Ratio 1: " << propValues[db_shearRatio  ] << "\n"
	<< "Ratio 2: " << propValues[db_shearRatio+1] << "\n"
	<< "Ratio 3: " << propValues[db_shearRatio+2] << "\n"
	<< "Total:   " << visFrac << "\n";
    throw std::runtime_error(msg.str());
  } // if

  // Loop over number of Maxwell models.
  for (int imodel =0; imodel < numMaxwellModels; ++imodel) {
    PylithScalar muRatio = dbValues[db_shearRatio + imodel];
    PylithScalar viscosity = dbValues[db_viscosity + imodel];
    PylithScalar muFac = muRatio*mu;
    PylithScalar maxwellTime = pylith::PYLITH_MAXSCALAR;
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
    propValues[p_shearRatio + imodel] = muRatio;
    propValues[p_maxwellTime + imodel] = maxwellTime;
  } // for

  PetscLogFlops(6 + 3 * numMaxwellModels);
} // _dbToProperties

// ----------------------------------------------------------------------
// Nondimensionalize properties.
void
pylith::materials::GenMaxwellIsotropic3D::_nondimProperties(PylithScalar* const values,
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
  values[p_muEff] = 
    _normalizer->nondimensionalize(values[p_muEff], pressureScale);
  values[p_lambdaEff] = 
    _normalizer->nondimensionalize(values[p_lambdaEff], pressureScale);
  const int numMaxwellModels = _GenMaxwellIsotropic3D::numMaxwellModels;
  _normalizer->nondimensionalize(&values[p_maxwellTime],
				 numMaxwellModels, timeScale);
  
  PetscLogFlops(3+1*numMaxwellModels);
} // _nondimProperties

// ----------------------------------------------------------------------
// Dimensionalize properties.
void
pylith::materials::GenMaxwellIsotropic3D::_dimProperties(PylithScalar* const values,
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
  values[p_muEff] = 
    _normalizer->dimensionalize(values[p_muEff], pressureScale);
  values[p_lambdaEff] = 
    _normalizer->dimensionalize(values[p_lambdaEff], pressureScale);
  const int numMaxwellModels = _GenMaxwellIsotropic3D::numMaxwellModels;
  _normalizer->dimensionalize(&values[p_maxwellTime],
			      numMaxwellModels, timeScale);

  PetscLogFlops(3+1*numMaxwellModels);
} // _dimProperties

// ----------------------------------------------------------------------
// Compute initial state variables from values in spatial database.
void
pylith::materials::GenMaxwellIsotropic3D::_dbToStateVars(
					PylithScalar* const stateValues,
					const scalar_array& dbValues)
{ // _dbToStateVars
  assert(0 != stateValues);
  const int numDBValues = dbValues.size();
  assert(_GenMaxwellIsotropic3D::numDBStateVars == numDBValues);

  const int numMaxwellModels = _GenMaxwellIsotropic3D::numMaxwellModels;
  const int totalSize = (1 + numMaxwellModels) * _tensorSize;
  assert(totalSize == _numVarsQuadPt);
  assert(totalSize == numDBValues);
  for (int i=0; i < totalSize; ++i)
    stateValues[i] = dbValues[i];

  PetscLogFlops(0);
} // _dbToStateVars

// ----------------------------------------------------------------------
// Compute density at location from properties.
void
pylith::materials::GenMaxwellIsotropic3D::_calcDensity(PylithScalar* const density,
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
pylith::materials::GenMaxwellIsotropic3D::_calcStressElastic(
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
  assert(_GenMaxwellIsotropic3D::tensorSize == stressSize);
  assert(0 != properties);
  assert(_numPropsQuadPt == numProperties);
  assert(0 != stateVars);
  assert(_numVarsQuadPt == numStateVars);
  assert(0 != totalStrain);
  assert(_GenMaxwellIsotropic3D::tensorSize == strainSize);
  assert(0 != initialStress);
  assert(_GenMaxwellIsotropic3D::tensorSize == initialStressSize);
  assert(0 != initialStrain);
  assert(_GenMaxwellIsotropic3D::tensorSize == initialStrainSize);

  const PylithScalar mu = properties[p_muEff];
  const PylithScalar lambda = properties[p_lambdaEff];
  const PylithScalar mu2 = 2.0 * mu;

  // :TODO: Need to consider initial state variables????
  const PylithScalar e11 = totalStrain[0] - initialStrain[0];
  const PylithScalar e22 = totalStrain[1] - initialStrain[1];
  const PylithScalar e33 = totalStrain[2] - initialStrain[2];
  const PylithScalar e12 = totalStrain[3] - initialStrain[3];
  const PylithScalar e23 = totalStrain[4] - initialStrain[4];
  const PylithScalar e13 = totalStrain[5] - initialStrain[5];
  
  const PylithScalar s123 = lambda * (e11 + e22 + e33);

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
pylith::materials::GenMaxwellIsotropic3D::_calcStressViscoelastic(
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
  assert(_GenMaxwellIsotropic3D::tensorSize == stressSize);
  assert(0 != properties);
  assert(_numPropsQuadPt == numProperties);
  assert(0 != stateVars);
  assert(_numVarsQuadPt == numStateVars);
  assert(0 != totalStrain);
  assert(_GenMaxwellIsotropic3D::tensorSize == strainSize);
  assert(0 != initialStress);
  assert(_GenMaxwellIsotropic3D::tensorSize == initialStressSize);
  assert(0 != initialStrain);
  assert(_GenMaxwellIsotropic3D::tensorSize == initialStrainSize);

  const int numMaxwellModels = _GenMaxwellIsotropic3D::numMaxwellModels;
  const int tensorSize = _tensorSize;

  const PylithScalar mu = properties[p_muEff];
  const PylithScalar lambda = properties[p_lambdaEff];
  const PylithScalar muRatio[numMaxwellModels] = {
    properties[p_shearRatio  ],
    properties[p_shearRatio+1],
    properties[p_shearRatio+2]
  };

  const PylithScalar mu2 = 2.0 * mu;
  const PylithScalar bulkModulus = lambda + mu2/3.0;

  // Initial stress and strain values
  const PylithScalar meanStrainInitial =
    (initialStrain[0] + initialStrain[1] + initialStrain[2])/3.0;
  const PylithScalar meanStressInitial =
    (initialStress[0] + initialStress[1] + initialStress[2])/3.0;
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
  // Mean stress and strain for time t + dt
  const PylithScalar meanStrainTpdt = (totalStrain[0] +
				 totalStrain[1] +
				 totalStrain[2]) / 3.0;
  
  const PylithScalar meanStressTpdt = 3.0 * bulkModulus *
    (meanStrainTpdt - meanStrainInitial) + meanStressInitial;

  const PylithScalar diag[] = { 1.0, 1.0, 1.0, 0.0, 0.0, 0.0 };

  PylithScalar visFrac = 0.0;
  for (int imodel=0; imodel < numMaxwellModels; ++imodel) 
    visFrac += muRatio[imodel];
  assert(visFrac <= 1.0);
  const PylithScalar elasFrac = 1.0 - visFrac;

  PetscLogFlops(23 + numMaxwellModels);

  // Get viscous strains
  if (computeStateVars) {
    _computeStateVars(stateVars, numStateVars,
		      properties, numProperties,
		      totalStrain, strainSize,
		      initialStress, initialStressSize,
		      initialStrain, initialStrainSize);
  } else {
    int index = 0;
    for (int iComp=0; iComp < tensorSize; ++iComp)
      _viscousStrain[index++] = stateVars[s_viscousStrain1+iComp];
    for (int iComp=0; iComp < tensorSize; ++iComp)
      _viscousStrain[index++] = stateVars[s_viscousStrain2+iComp];
    for (int iComp=0; iComp < tensorSize; ++iComp)
      _viscousStrain[index++] = stateVars[s_viscousStrain3+iComp];
  } // else

  // Compute new stresses
  PylithScalar devStrainTpdt = 0.0;
  PylithScalar devStressTpdt = 0.0;
  for (int iComp=0; iComp < tensorSize; ++iComp) {
    devStrainTpdt = totalStrain[iComp] - diag[iComp] * meanStrainTpdt -
      devStrainInitial[iComp];
    devStressTpdt = elasFrac * devStrainTpdt;
    for (int model=0; model < numMaxwellModels; ++model) {
      devStressTpdt += muRatio[model] *
	_viscousStrain[model * tensorSize+iComp];
    } // for

    devStressTpdt = mu2 * devStressTpdt + devStressInitial[iComp];
    stress[iComp] = diag[iComp] * meanStressTpdt + devStressTpdt;
  } // for

  PetscLogFlops((9 + 3 * numMaxwellModels) * tensorSize);
} // _calcStressViscoelastic

// ----------------------------------------------------------------------
// Compute derivative of elasticity matrix at location from properties.
void
pylith::materials::GenMaxwellIsotropic3D::_calcElasticConstsElastic(
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
  assert(_GenMaxwellIsotropic3D::numElasticConsts == numElasticConsts);
  assert(0 != properties);
  assert(_numPropsQuadPt == numProperties);
  assert(0 != stateVars);
  assert(_numVarsQuadPt == numStateVars);
  assert(0 != totalStrain);
  assert(_GenMaxwellIsotropic3D::tensorSize == strainSize);
  assert(0 != initialStress);
  assert(_GenMaxwellIsotropic3D::tensorSize == initialStressSize);
  assert(0 != initialStrain);
  assert(_GenMaxwellIsotropic3D::tensorSize == initialStrainSize);
 
  const PylithScalar mu = properties[p_muEff];
  const PylithScalar lambda = properties[p_lambdaEff];

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
// as a viscoelastic material.
void
pylith::materials::GenMaxwellIsotropic3D::_calcElasticConstsViscoelastic(
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
  assert(_GenMaxwellIsotropic3D::numElasticConsts == numElasticConsts);
  assert(0 != properties);
  assert(_numPropsQuadPt == numProperties);
  assert(0 != stateVars);
  assert(_numVarsQuadPt == numStateVars);
  assert(0 != totalStrain);
  assert(_GenMaxwellIsotropic3D::tensorSize == strainSize);
  assert(0 != initialStress);
  assert(_GenMaxwellIsotropic3D::tensorSize == initialStressSize);
  assert(0 != initialStrain);
  assert(_GenMaxwellIsotropic3D::tensorSize == initialStrainSize);

  const int numMaxwellModels = _GenMaxwellIsotropic3D::numMaxwellModels;

  const PylithScalar mu = properties[p_muEff];
  const PylithScalar lambda = properties[p_lambdaEff];
  const PylithScalar mu2 = 2.0 * mu;
  const PylithScalar bulkModulus = lambda + mu2 / 3.0;

  // Compute viscous contribution.
  PylithScalar visFac = 0.0;
  PylithScalar visFrac = 0.0;
  PylithScalar shearRatio = 0.0;
  for (int imodel = 0; imodel < numMaxwellModels; ++imodel) {
    shearRatio = properties[p_shearRatio + imodel];
    PylithScalar maxwellTime = pylith::PYLITH_MAXSCALAR;
    visFrac += shearRatio;
    if (shearRatio != 0.0) {
      maxwellTime = properties[p_maxwellTime + imodel];
      visFac +=
	shearRatio*ViscoelasticMaxwell::viscousStrainParam(_dt, maxwellTime);
    } // if
  } // for
  PylithScalar elasFrac = 1.0 - visFrac;
  PylithScalar shearFac = elasFrac + visFac;

  elasticConsts[ 0] = bulkModulus + 4.0*mu/3.0 * shearFac; // C1111
  elasticConsts[ 1] = bulkModulus - 2.0*mu/3.0 * shearFac; // C1122
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
  elasticConsts[21] = 2.0*mu*shearFac; // C1212
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
// Update state variables.
void
pylith::materials::GenMaxwellIsotropic3D::_updateStateVarsElastic(
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
  assert(_GenMaxwellIsotropic3D::tensorSize == strainSize);
  assert(0 != initialStress);
  assert(_GenMaxwellIsotropic3D::tensorSize == initialStressSize);
  assert(0 != initialStrain);
  assert(_GenMaxwellIsotropic3D::tensorSize == initialStrainSize);

  const int tensorSize = _tensorSize;

  const PylithScalar strainTpdt[] = {totalStrain[0] - initialStrain[0],
			       totalStrain[1] - initialStrain[1],
			       totalStrain[2] - initialStrain[2],
			       totalStrain[3] - initialStrain[3],
			       totalStrain[4] - initialStrain[4],
			       totalStrain[5] - initialStrain[5]};
  const PylithScalar meanStrainTpdt =
    (strainTpdt[0] + strainTpdt[1] + strainTpdt[2])/3.0;

  const PylithScalar diag[] = { 1.0, 1.0, 1.0, 0.0, 0.0, 0.0 };

  // Update total strain
  for (int iComp=0; iComp < tensorSize; ++iComp)
    stateVars[s_totalStrain+iComp] = totalStrain[iComp];

  // Initialize all viscous strains to deviatoric elastic strains.
  PylithScalar devStrain = 0.0;
  for (int iComp=0; iComp < tensorSize; ++iComp) {
    devStrain = strainTpdt[iComp] - diag[iComp] * meanStrainTpdt;
    // Maxwell model 1
    stateVars[s_viscousStrain1+iComp] = devStrain;
    // Maxwell model 2
    stateVars[s_viscousStrain2+iComp] = devStrain;
    // Maxwell model 3
    stateVars[s_viscousStrain3+iComp] = devStrain;
  } // for
  PetscLogFlops(9 + 2 * tensorSize);

  _needNewJacobian = true;
} // _updateStateVarsElastic

// ----------------------------------------------------------------------
// Update state variables.
void
pylith::materials::GenMaxwellIsotropic3D::_updateStateVarsViscoelastic(
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
  assert(_GenMaxwellIsotropic3D::tensorSize == strainSize);
  assert(0 != initialStress);
  assert(_GenMaxwellIsotropic3D::tensorSize == initialStressSize);
  assert(0 != initialStrain);
  assert(_GenMaxwellIsotropic3D::tensorSize == initialStrainSize);

  _computeStateVars(stateVars, numStateVars,
		    properties, numProperties,
		    totalStrain, strainSize,
		    initialStress, initialStressSize,
		    initialStrain, initialStrainSize);

  const int tensorSize = _tensorSize;

  // Total strain
  for (int iComp=0; iComp < tensorSize; ++iComp)
    stateVars[s_totalStrain+iComp] = totalStrain[iComp];

  // Viscous strain
  int index = 0;
  for (int iComp=0; iComp < tensorSize; ++iComp)
    stateVars[s_viscousStrain1+iComp] = _viscousStrain[index++];
  for (int iComp=0; iComp < tensorSize; ++iComp)
    stateVars[s_viscousStrain2+iComp] = _viscousStrain[index++];
  for (int iComp=0; iComp < tensorSize; ++iComp)
    stateVars[s_viscousStrain3+iComp] = _viscousStrain[index++];

  _needNewJacobian = false;
} // _updateStateVarsViscoelastic

// ----------------------------------------------------------------------
// Get stable time step for implicit time integration.
PylithScalar
pylith::materials::GenMaxwellIsotropic3D::_stableTimeStepImplicit(
					   const PylithScalar* properties,
					   const int numProperties,
					   const PylithScalar* stateVars,
					   const int numStateVars) const
{ // _stableTimeStepImplicit
  assert(0 != properties);
  assert(_numPropsQuadPt == numProperties);
  assert(0 != stateVars);
  assert(_numVarsQuadPt == numStateVars);

  PylithScalar dtStable = pylith::PYLITH_MAXSCALAR;

  const int numMaxwellModels = _GenMaxwellIsotropic3D::numMaxwellModels;
  for (int i=0; i < numMaxwellModels; ++i) {
    const PylithScalar maxwellTime = properties[p_maxwellTime+i];
    const PylithScalar dt = 0.2*maxwellTime;
    if (dt < dtStable)
      dtStable = dt;
  } // for

  return dtStable;
} // _stableTimeStepImplicit

// ----------------------------------------------------------------------
// Get stable time step for explicit time integration.
PylithScalar
pylith::materials::GenMaxwellIsotropic3D::_stableTimeStepExplicit(const PylithScalar* properties,
								  const int numProperties,
								  const PylithScalar* stateVars,
								  const int numStateVars,
								  const double minCellWidth) const
{ // _stableTimeStepExplicit
  assert(properties);
  assert(_numPropsQuadPt == numProperties);
 
  const PylithScalar mu = properties[p_muEff];
  const PylithScalar lambda = properties[p_lambdaEff];
  const PylithScalar density = properties[p_density];

  assert(density > 0.0);
  const PylithScalar vp = sqrt((lambda + 2*mu) / density);

  const PylithScalar dtStable = minCellWidth / vp;
  return dtStable;
} // _stableTimeStepExplicit


// ----------------------------------------------------------------------
// Compute viscous strain for current time step.
void
pylith::materials::GenMaxwellIsotropic3D::_computeStateVars(
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
  assert(_GenMaxwellIsotropic3D::tensorSize == strainSize);
  assert(0 != initialStress);
  assert(_GenMaxwellIsotropic3D::tensorSize == initialStressSize);
  assert(0 != initialStrain);
  assert(_GenMaxwellIsotropic3D::tensorSize == initialStrainSize);

  const int tensorSize = _tensorSize;
  const int numMaxwellModels = _GenMaxwellIsotropic3D::numMaxwellModels;

  const PylithScalar muRatio[numMaxwellModels] = {
    properties[p_shearRatio  ],
    properties[p_shearRatio+1],
    properties[p_shearRatio+2]
  };
  const PylithScalar maxwellTime[numMaxwellModels] = {
    properties[p_maxwellTime  ],
    properties[p_maxwellTime+1],
    properties[p_maxwellTime+2]
  };

  // :TODO: Need to account for initial values for state variables
  const PylithScalar meanStrainTpdt =
    (totalStrain[0] + totalStrain[1] + totalStrain[2])/3.0;
  const PylithScalar diag[] = { 1.0, 1.0, 1.0, 0.0, 0.0, 0.0 };

  const PylithScalar meanStrainT = 
    ( stateVars[s_totalStrain+0] +
      stateVars[s_totalStrain+1] +
      stateVars[s_totalStrain+2] ) / 3.0;
  
  PetscLogFlops(6);

  // Compute Prony series terms
  scalar_array dq(numMaxwellModels);
  dq = 0.0;
  for (int i=0; i < numMaxwellModels; ++i)
    if (muRatio[i] != 0.0)
      dq[i] = ViscoelasticMaxwell::viscousStrainParam(_dt, maxwellTime[i]);

  // Compute new viscous strains
  PylithScalar devStrainTpdt = 0.0;
  PylithScalar devStrainT = 0.0;
  PylithScalar deltaStrain = 0.0;
  
  for (int iComp=0; iComp < tensorSize; ++iComp) {
    devStrainTpdt = totalStrain[iComp] - diag[iComp] * meanStrainTpdt;
    devStrainT = stateVars[s_totalStrain+iComp] - diag[iComp] * meanStrainT;
    deltaStrain = devStrainTpdt - devStrainT;

    // Maxwell model 1
    int imodel = 0;
    if (0.0 != muRatio[imodel]) {
      _viscousStrain[imodel * tensorSize+iComp] = 
	exp(-_dt/maxwellTime[imodel]) *
	stateVars[s_viscousStrain1 + iComp] + dq[imodel] * deltaStrain;
      PetscLogFlops(6);
    } // if

    // Maxwell model 2
    imodel = 1;
    if (0.0 != muRatio[imodel]) {
      _viscousStrain[imodel*tensorSize+iComp] =
	exp(-_dt/maxwellTime[imodel]) *
	stateVars[s_viscousStrain2 + iComp] + dq[imodel] * deltaStrain;
      PetscLogFlops(6);
    } // if

    // Maxwell model 3
    imodel = 2;
    if (0.0 != muRatio[imodel]) {
      _viscousStrain[imodel*tensorSize+iComp] =
	exp(-_dt/maxwellTime[imodel]) *
	stateVars[s_viscousStrain3 + iComp] + dq[imodel] * deltaStrain;
      PetscLogFlops(6);
    } // if

  } // for

  PetscLogFlops(5 * tensorSize);
} // _computeStateVars


// End of file 
