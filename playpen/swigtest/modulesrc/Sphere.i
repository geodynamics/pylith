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

namespace geometry {

class Sphere : public Shape
{ // Sphere
public :
  
  Sphere(void);
  ~Sphere(void);

  void radius(const double value);
  double radius(void) const;

  void view(void) const;

}; // Sphere

} // namespace geometry


// End of file
