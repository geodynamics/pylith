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
#include "Sphere.hh"
#include "Scene.hh"

int
main(int argc,
     char** argv)
{
  geometry::Bar bar;
  bar.length(4.0);
  bar.width(2.0);
  bar.color("yellow");
  bar.id(0);

  geometry::Sphere sphere;
  sphere.radius(3.0);
  sphere.color("green");
  sphere.id(1);

  geometry::Scene scene;
  scene.bar(&bar);
  scene.sphere(&sphere);
  scene.view();

  return 0;
}


// End of file
