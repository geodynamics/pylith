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

#if !defined(scene_hh)
#define scene_hh

namespace geometry {
  class Scene;

  class Bar;
  class Sphere;
} // geometry

class geometry::Scene
{ // Scene
public :
  
  Scene(void);
  ~Scene(void);

  void bar(Bar* const value);
  const Bar* bar(void) const;

  void sphere(Sphere* const value);
  const Sphere* sphere(void) const;

  void view(void) const;

private :

  Bar* _bar;
  Sphere* _sphere;

}; // Scene

#endif // scene_hh


// End of file
