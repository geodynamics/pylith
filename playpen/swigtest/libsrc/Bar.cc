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


#include "Bar.hh"

#include <iostream>
#include <stdexcept>

geometry::Bar::Bar(void) :
  _length(0),
  _width(0)
{}
geometry::Bar::~Bar(void)
{}

void
geometry::Bar::length(const double value) {
  if (value <= 0)
    throw std::runtime_error("Length must be positive.");
  _length = value;
}
double
geometry::Bar::length(void) const {
  return _length;
}

void
geometry::Bar::width(const double value) {
  if (value <= 0)
    throw std::runtime_error("Width must be positive.");
  _width = value;
}
double
geometry::Bar::width(void) const {
  return _width;
}

void
geometry::Bar::view(void) const {
  std::cout << "Bar:\n"
	    << "  length: " << _length << "\n"
	    << "  width: " << _width << "\n"
	    << "  color: " << _color << "\n"
	    << "  id: " << _id << "\n"
	    << std::endl;
}


// End of file
