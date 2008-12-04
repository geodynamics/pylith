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

#include <portinfo>

#include "Field.hh" // implementation of class methods

#include "pylith/utils/array.hh" // USES double_array

#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include <stdexcept> // USES std::runtime_error
#include <sstream> // USES std::ostringstream
#include <cassert> // USES assert()

// ----------------------------------------------------------------------
// Default constructor.
pylith::topology::Field::Field(const ALE::Obj<SieveMesh>& mesh) :
  _scale(1.0),
  _name("unknown"),
  _spaceDim(0),
  _vecFieldType(OTHER),
  _dimensionsOkay(false)
{ // constructor
  assert(!mesh.isNull());

  _section = new SieveMesh::real_section_type(mesh->comm(), mesh->debug());
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::topology::Field::~Field(void)
{ // destructor
} // destructor

// ----------------------------------------------------------------------
// Dimensionalize field.
void
pylith::topology::Field::dimensionalize(void)
{ // dimensionalize
  if (!_dimensionsOkay) {
    std::ostringstream msg;
    msg << "Cannot dimensionalize field '" << _name << "' because the flag "
	<< "has been set to keep field nondimensional.";
    throw std::runtime_error(msg.str());
  } // if

  typedef SieveMesh::real_section_type::chart_type chart_type;

  assert(!_section.isNull());
  const chart_type& chart = _section->getChart();
  const chart_type::const_iterator chartEnd = chart.end();

  // Assume fiber dimension is uniform
  const int fiberDim = _section->getFiberDimension(*chart.begin());
  double_array values(fiberDim);

  spatialdata::units::Nondimensional normalizer;

  for (chart_type::const_iterator c_iter = chart.begin();
       c_iter != chartEnd;
       ++c_iter) {
    assert(fiberDim == _section->getFiberDimension(*c_iter));

    _section->restrictPoint(*c_iter, &values[0], values.size());
    normalizer.dimensionalize(&values[0], values.size(), _scale);
    _section->updatePoint(*c_iter, &values[0]);
  } // for
  
} // dimensionalize


// End of file 
