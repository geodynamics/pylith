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


#include "Scene.hh"

#include "Bar.hh"
#include "Sphere.hh"

#include <iostream>

geometry::Scene::Scene(void) :
  _bar(0),
  _sphere(0)
{}
geometry::Scene::~Scene(void)
{}

void
geometry::Scene::bar(Bar* const value) {
  _bar = value;
}
const geometry::Bar*
geometry::Scene::bar(void) const {
  return _bar;
}

void
geometry::Scene::sphere(Sphere* const value) {
  _sphere = value;
}
const geometry::Sphere*
geometry::Scene::sphere(void) const {
  return _sphere;
}

void
geometry::Scene::view(void) const {
  std::cout << "Scene:\n";
  if (0 != _bar)
    _bar->view();
  else
    std::cout << "  Bar: None\n";
  if (0 != _sphere)
    _sphere->view();
  else
    std::cout << "  Sphere: None\n";
}


// End of file
