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
  density(0),
  strain(0),
  stress(0),
  elasticConsts(0)
{ // constructor
} // constructor

pylith::materials::ElasticMaterialData::~ElasticMaterialData(void)
{}


// End of file
