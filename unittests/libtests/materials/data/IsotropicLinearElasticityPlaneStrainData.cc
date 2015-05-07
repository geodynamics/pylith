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

// ----------------------------------------------------------------------
// Constructor
pylith::materials::IsotropicLinearElasticityPlaneStrainData::IsotropicLinearElasticityPlaneStrainData(void) :
  filenameMesh(0),
  useInertia(false),
  useBodyForce(false),
  numSolnFields(0),
  discretizations(NULL),
  numAuxFields(0),
  auxSubfields(NULL),
  auxFields(NULL),
  lengthScale(1.0e+3),
  timeScale(2.0),
  pressureScale(2.25e+10),
  densityScale(0)
{ // constructor
  const PylithScalar velScale = lengthScale / timeScale;
  densityScale = pressureScale / (velScale*velScale);
} // constructor


// ----------------------------------------------------------------------
// Destructor
pylith::materials::IsotropicLinearElasticityPlaneStrainData::~IsotropicLinearElasticityPlaneStrainData(void)
{ // destructor
} // destructor

// End of file
