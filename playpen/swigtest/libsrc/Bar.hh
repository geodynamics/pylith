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

#if !defined(bar_hh)
#define bar_hh

#include "Shape.hh"

namespace geometry {
  class Bar;
} // geometry

class geometry::Bar : public Shape
{ // Bar
public :
  
  Bar(void);
  ~Bar(void);

  void length(const double value);
  double length(void) const;

  void width(const double value);
  double width(void) const;

  void view(void) const;

private :

  double _length;
  double _width;

}; // Bar

#endif // bar_hh


// End of file
