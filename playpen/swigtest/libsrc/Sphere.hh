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

#if !defined(sphere_hh)
#define sphere_hh

#include "Shape.hh"

namespace geometry {
  class Sphere;
} // geometry

class geometry::Sphere : public Shape
{ // Sphere
public :
  
  Sphere(void);
  ~Sphere(void);

  void radius(const double value);
  double radius(void) const;

  void view(void) const;

private :

  double _radius;

}; // Sphere

#endif // sphere_hh


// End of file
