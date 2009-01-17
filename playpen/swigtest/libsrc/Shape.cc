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


#include "Shape.hh"

#include <stdexcept>

geometry::Shape::Shape(void) :
  _id(0)
{}
geometry::Shape::~Shape(void)
{}

void
geometry::Shape::color(const char* value) {
  _color = value;
}
const char*
geometry::Shape::color(void) const {
  return _color.c_str();
}

void
geometry::Shape::id(const int value) {
  _id = value;
}
int
geometry::Shape::id(void) const {
  return _id;
}


// End of file
