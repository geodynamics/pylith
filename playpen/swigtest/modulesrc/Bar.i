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

class Bar : public Shape
{ // Bar
public :
  
  Bar(void);
  ~Bar(void);

  void length(const double value);
  double length(void) const;

  void width(const double value);
  double width(void) const;

  void view(void) const;

}; // Bar


} // namespace geometry

// End of file
