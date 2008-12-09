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

class Shape
{ // Shape
public :
  
  Shape(void);
  virtual ~Shape(void);

  void color(const char* value);
  const char* color(void) const;

  void id(const int value);
  int id(void) const;

  virtual void view(void) const = 0;

}; // Shape

} // namespace geometry

// End of file
