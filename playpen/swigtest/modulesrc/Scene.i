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

}; // Scene

} // namespace geometry


// End of file
