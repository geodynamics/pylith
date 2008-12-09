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

#if !defined(shape_hh)
#define shape_hh

#include <string>

namespace geometry {
  class Shape;
} // geometry

class geometry::Shape
{ // Shape
public :
  
  Shape(void);
  virtual ~Shape(void);

  void color(const char* value);
  const char* color(void) const;

  void id(const int value);
  int id(void) const;

  virtual void view(void) const;

protected :

  std::string _color;
  int _id;

}; // Shape

#endif // shape_hh


// End of file
