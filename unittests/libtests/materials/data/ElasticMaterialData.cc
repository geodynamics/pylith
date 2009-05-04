// -*- C++ -*-
//
// ======================================================================
//
//                           Brad T. Aagaard
//                        U.S. Geological Survey
//
// {LicenseText}
//
// ======================================================================
//

#include "ElasticMaterialData.hh"

pylith::materials::ElasticMaterialData::ElasticMaterialData(void) :
  dtStableImplicit(0.0),
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
