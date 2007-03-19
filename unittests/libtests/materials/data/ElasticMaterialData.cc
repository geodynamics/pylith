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
  numLocs(0),
  density(0),
  elasticConsts(0)
{ // constructor
} // constructor

pylith::materials::ElasticMaterialData::~ElasticMaterialData(void)
{}


// End of file
