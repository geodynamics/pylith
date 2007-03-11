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

#include "ElasticMaterial3D.hh" // implementation of object methods

#include <assert.h> // USES assert()

// ----------------------------------------------------------------------
const int pylith::materials::ElasticMaterial3D::NUMELASTCONSTS = 21;

// ----------------------------------------------------------------------
// Default constructor.
pylith::materials::ElasticMaterial3D::ElasticMaterial3D(void)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::materials::ElasticMaterial3D::~ElasticMaterial3D(void)
{ // destructor
} // destructor

// ----------------------------------------------------------------------
// Copy constructor.
pylith::materials::ElasticMaterial3D::ElasticMaterial3D(const ElasticMaterial3D& m) :
  Material(m)
{ // copy constructor
} // copy constructor

// ----------------------------------------------------------------------
// Get section with inertia values.
const ALE::Obj<ALE::Mesh::real_section_type>&
pylith::materials::ElasticMaterial3D::density(void)
{ // density
} // density

// ----------------------------------------------------------------------
// Get section with elasticity constants.
const ALE::Obj<ALE::Mesh::real_section_type>&
pylith::materials::ElasticMaterial3D::elasticityConstants(void)
{ // elasticityConstants
} // elasticityConstants


// End of file 
