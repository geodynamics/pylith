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

#include "ElasticIsotropic3D.hh" // implementation of object methods

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
    namespace _ElasticIsotropic3D {

      // Dimension of material.
      const int dimension = 3;

      // Number of entries in stress tensor.
      const int tensorSize = 6;

      // Number of elastic constants (for general 3-D elastic material)
      const int numElasticConsts = 36;

      // Number of physical properties.
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
      
    } // _ElasticIsotropic3D
  } // materials
} // pylith

// Indices of physical properties
const int pylith::materials::ElasticIsotropic3D::p_density = 0;

const int pylith::materials::ElasticIsotropic3D::p_mu = 
  pylith::materials::ElasticIsotropic3D::p_density + 1;

const int pylith::materials::ElasticIsotropic3D::p_lambda = 
  pylith::materials::ElasticIsotropic3D::p_mu + 1;

// Indices of database values (order must match dbProperties)
const int pylith::materials::ElasticIsotropic3D::db_density = 0;

const int pylith::materials::ElasticIsotropic3D::db_vs = 
  pylith::materials::ElasticIsotropic3D::db_density + 1;

const int pylith::materials::ElasticIsotropic3D::db_vp = 
  pylith::materials::ElasticIsotropic3D::db_vs + 1;

// ----------------------------------------------------------------------
// Default constructor.
pylith::materials::ElasticIsotropic3D::ElasticIsotropic3D(void) :
  ElasticMaterial(_ElasticIsotropic3D::dimension,
		  _ElasticIsotropic3D::tensorSize,
		  _ElasticIsotropic3D::numElasticConsts,
		  Metadata(_ElasticIsotropic3D::properties,
			   _ElasticIsotropic3D::numProperties,
			   _ElasticIsotropic3D::dbProperties,
			   _ElasticIsotropic3D::numDBProperties,
			   0, 0,
			   0, 0))
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::materials::ElasticIsotropic3D::~ElasticIsotropic3D(void)
{ // destructor
} // destructor

// ----------------------------------------------------------------------
// Compute properties from values in spatial database.
void
pylith::materials::ElasticIsotropic3D::_dbToProperties(
					   PylithScalar* const propValues,
					   const scalar_array& dbValues)
{ // _dbToProperties
  assert(0 != propValues);
  const int numDBValues = dbValues.size();
  assert(_ElasticIsotropic3D::numDBProperties == numDBValues);

  const PylithScalar density = dbValues[db_density];
  const PylithScalar vs = dbValues[db_vs];
  const PylithScalar vp = dbValues[db_vp];
 
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

  propValues[p_density] = density;
  propValues[p_mu] = mu;
  propValues[p_lambda] = lambda;

  PetscLogFlops(6);
} // _dbToProperties

// ----------------------------------------------------------------------
// Nondimensionalize properties.
void
pylith::materials::ElasticIsotropic3D::_nondimProperties(PylithScalar* const values,
							 const int nvalues) const
{ // _nondimProperties
  assert(0 != _normalizer);
  assert(0 != values);
  assert(nvalues == _ElasticIsotropic3D::numProperties);

  const PylithScalar densityScale = _normalizer->densityScale();
  const PylithScalar pressureScale = _normalizer->pressureScale();

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
pylith::materials::ElasticIsotropic3D::_dimProperties(PylithScalar* const values,
						      const int nvalues) const
{ // _dimProperties
  assert(0 != _normalizer);
  assert(0 != values);
  assert(nvalues == _ElasticIsotropic3D::numProperties);

  const PylithScalar densityScale = _normalizer->densityScale();
  const PylithScalar pressureScale = _normalizer->pressureScale();

  values[p_density] = 
    _normalizer->dimensionalize(values[p_density], densityScale);
  values[p_mu] = 
    _normalizer->dimensionalize(values[p_mu], pressureScale);
  values[p_lambda] = 
    _normalizer->dimensionalize(values[p_lambda], pressureScale);

  PetscLogFlops(3);
} // _dimProperties

// ----------------------------------------------------------------------
// Compute stress tensor at location from properties.
void
pylith::materials::ElasticIsotropic3D::_calcStress(PylithScalar* const stress,
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
{ // _calcStress
  assert(stress);
  assert(_ElasticIsotropic3D::tensorSize == stressSize);
  assert(properties);
  assert(_numPropsQuadPt == numProperties);
  assert(0 == numStateVars);
  assert(totalStrain);
  assert(_ElasticIsotropic3D::tensorSize == strainSize);
  assert(initialStress);
  assert(_ElasticIsotropic3D::tensorSize == initialStressSize);
  assert(initialStrain);
  assert(_ElasticIsotropic3D::tensorSize == initialStrainSize);

  const PylithScalar mu = properties[p_mu];
  const PylithScalar lambda = properties[p_lambda];

  const PylithScalar mu2 = 2.0*mu;

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
} // _calcStress

// ----------------------------------------------------------------------
// Compute derivative of elasticity matrix at location from properties.
void
pylith::materials::ElasticIsotropic3D::_calcElasticConsts(
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
{ // _calcElasticConsts
  assert(0 != elasticConsts);
  assert(_ElasticIsotropic3D::numElasticConsts == numElasticConsts);
  assert(0 != properties);
  assert(_numPropsQuadPt == numProperties);
  assert(0 == numStateVars);
  assert(0 != totalStrain);
  assert(_ElasticIsotropic3D::tensorSize == strainSize);
  assert(0 != initialStress);
  assert(_ElasticIsotropic3D::tensorSize == initialStressSize);
  assert(0 != initialStrain);
  assert(_ElasticIsotropic3D::tensorSize == initialStrainSize);
 
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
} // _calcElasticConsts

// ----------------------------------------------------------------------
// Get stable time step for implicit time integration.
PylithScalar
pylith::materials::ElasticIsotropic3D::stableTimeStepImplicit(const topology::Mesh& mesh,
							      topology::Field* field) {
  return ElasticMaterial::_stableTimeStepImplicitMax(mesh, field);
} // stableTimeStepImplicitMax

// ----------------------------------------------------------------------
// Get stable time step for implicit time integration.
PylithScalar
pylith::materials::ElasticIsotropic3D::_stableTimeStepImplicit(
				     const PylithScalar* properties,
				     const int numProperties,
				     const PylithScalar* stateVars,
				     const int numStateVars) const
{ // _stableTimeStepImplicit
  return pylith::PYLITH_MAXSCALAR;
} // _stableTimeStepImplicit


// ----------------------------------------------------------------------
// Get stable time step for explicit time integration.
PylithScalar
pylith::materials::ElasticIsotropic3D::_stableTimeStepExplicit(const PylithScalar* properties,
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


// End of file 
