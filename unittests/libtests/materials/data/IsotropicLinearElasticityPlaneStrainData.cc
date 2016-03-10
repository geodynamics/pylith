// -*- C++ -*-
//
// ======================================================================
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
// ======================================================================
//

#include "IsotropicLinearElasticityPlaneStrainData.hh"


const int pylith::materials::IsotropicLinearElasticityPlaneStrainData::numKernelsResidual = 2;
const int pylith::materials::IsotropicLinearElasticityPlaneStrainData::numKernelsJacobian = 4;

// ----------------------------------------------------------------------
// Constructor
pylith::materials::IsotropicLinearElasticityPlaneStrainData::IsotropicLinearElasticityPlaneStrainData(void) :
  meshFilename(0),
  materialLabel(NULL),
  materialId(0),
  boundaryLabel(NULL),

  lengthScale(0),
  timeScale(0),
  pressureScale(0),
  densityScale(0),

  t(0.0),
  dt(0.0),
  
  solnDiscretizations(NULL),

  numAuxFields(0),
  auxFields(NULL),
  auxDiscretizations(NULL),
  filenameAuxFieldsDB(NULL),

  dimension(0),
  numSolnFields(0),

  kernelsRHSResidual(NULL),
  kernelsRHSJacobian(NULL),
  kernelsLHSResidual(NULL),
  kernelsLHSJacobianImplicit(NULL),
  kernelsLHSJacobianExplicit(NULL),

  useInertia(false),
  useBodyForce(false)

{ // constructor
} // constructor


// ----------------------------------------------------------------------
// Destructor
pylith::materials::IsotropicLinearElasticityPlaneStrainData::~IsotropicLinearElasticityPlaneStrainData(void)
{ // destructor
} // destructor

// End of file
