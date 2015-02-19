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

#include "ElasticMaterialData.hh"

pylith::materials::ElasticMaterialData::ElasticMaterialData(void) :
  dtStableImplicit(0.0),
  dtStableExplicit(0.0),
  density(0),
  strain(0),
  stress(0),
  elasticConsts(0),
  initialStress(0),
  initialStrain(0),
  stateVarsUpdated(0)
{ // constructor
} // constructor

pylith::materials::ElasticMaterialData::~ElasticMaterialData(void)
{}


// End of file
