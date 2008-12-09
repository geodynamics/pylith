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


#include "Sphere.hh"

#include <iostream>
#include <stdexcept>

geometry::Sphere::Sphere(void) :
  _radius(0)
{}
geometry::Sphere::~Sphere(void)
{}

void
geometry::Sphere::radius(const double value) {
  if (value <= 0)
    throw std::runtime_error("Radius must be positive.");
  _radius = value;
}
double
geometry::Sphere::radius(void) const {
  return _radius;
}

void
geometry::Sphere::view(void) const {
  std::cout << "Sphere:\n"
	    << "  radius: " << _radius << "\n"
	    << "  color: " << _color << "\n"
	    << "  id: " << _id << "\n"
	    << std::endl;
}


// End of file
