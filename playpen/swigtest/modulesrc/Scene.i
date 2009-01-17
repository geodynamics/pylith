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

class Bar;
class Sphere;

class Scene
{ // Scene
public :
  
  Scene(void);
  ~Scene(void);

  void bar(geometry::Bar* const value);
  const geometry::Bar* bar(void) const;

  void sphere(geometry::Sphere* const value);
  const geometry::Sphere* sphere(void) const;

  void view(void) const;

  %apply(double* IN_ARRAY1, int DIM1) {(const double* values, const int size)};
  void printData(const double* values,
		 const int size) const;  
  %clear (const double* values, const int size);

}; // Scene

} // namespace geometry


// End of file
