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
  _mesh(mesh),
  _scale(1.0),
  _name("unknown"),
  _spaceDim(0),
  _vecFieldType(OTHER),
  _dimensionsOkay(false)
{ // constructor
  assert(!mesh.isNull());

  _section = new SieveRealSection(mesh->comm(), mesh->debug());
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::topology::Field::~Field(void)
{ // destructor
} // destructor

// ----------------------------------------------------------------------
// Create section given atlas.
void
pylith::topology::Field::copyLayout(const Field& src)
{ // createSection
  _spaceDim = src._spaceDim;
  _vecFieldType = src._vecFieldType;

  const ALE::Obj<SieveRealSection>& srcSection = src.section();
  assert(!_section.isNull());

  _section->setAtlas(srcSection->getAtlas());
  _section->allocateStorage();
  _section->setBC(srcSection->getBC());
} // createSection

// ----------------------------------------------------------------------
// Clear variables associated with section.
void
pylith::topology::Field::clear(void)
{ // clear
  assert(!_section.isNull());
  _section->clear();

  _scale = 1.0;
  _vecFieldType = OTHER;
  _dimensionsOkay = false;
} // clear

// ----------------------------------------------------------------------
// Zero section values.
void
pylith::topology::Field::zero(void)
{ // zero
  assert(!_section.isNull());
  _section->zero();
} // zero

// ----------------------------------------------------------------------
// Complete section by assembling across processors.
void
pylith::topology::Field::complete(void)
{ // complete
  assert(!_section.isNull());
  ALE::Completion::completeSectionAdd(_mesh->getSendOverlap(),
				      _mesh->getRecvOverlap(), 
				      _section, _section);
} // complete

// ----------------------------------------------------------------------
// Copy field values and metadata.
void
pylith::topology::Field::copy(const Field& field)
{ // copy
  // Check compatibility of sections
  assert(!_section.isNull());
  assert(!field._section.isNull());
  const int srcSize = field._section->size();
  const int dstSize = _section->size();
  if (field._spaceDim != _spaceDim ||
      field._vecFieldType != _vecFieldType ||
      field._scale != _scale ||
      srcSize != dstSize) {
    std::ostringstream msg;

    msg << "Cannot copy values from section '" << field._name 
	<< "' to section '" << _name << "'. Sections are incompatible.\n"
	<< "  Source section:\n"
	<< "    space dim: " << field._spaceDim << "\n"
	<< "    vector field type: " << field._vecFieldType << "\n"
	<< "    scale: " << field._scale << "\n"
	<< "    size: " << srcSize
	<< "  Destination section:\n"
	<< "    space dim: " << _spaceDim << "\n"
	<< "    vector field type: " << _vecFieldType << "\n"
	<< "    scale: " << _scale << "\n"
	<< "    size: " << dstSize;
    throw std::runtime_error(msg.str());
  } // if

  // Copy values from field
  const SieveRealSection::chart_type& chart = _section->getChart();
  const SieveRealSection::chart_type::const_iterator chartEnd = chart.end();

  for (SieveRealSection::chart_type::const_iterator c_iter = chart.begin();
       c_iter != chartEnd;
       ++c_iter) {
    assert(field._section->getFiberDimension(*c_iter) ==
	   _section->getFiberDimension(*c_iter));
    _section->updatePoint(*c_iter, field._section->restrictPoint(*c_iter));
  } // for
} // copy

// ----------------------------------------------------------------------
// Add two fields, storing the result in one of the fields.
void
pylith::topology::Field::operator+=(const Field& field)
{ // operator+=
  // Check compatibility of sections
  assert(!_section.isNull());
  assert(!field._section.isNull());
  const int srcSize = field._section->size();
  const int dstSize = _section->size();
  if (field._spaceDim != _spaceDim ||
      field._vecFieldType != _vecFieldType ||
      field._scale != _scale ||
      srcSize != dstSize) {
    std::ostringstream msg;

    msg << "Cannot add values from section '" << field._name 
	<< "' to section '" << _name << "'. Sections are incompatible.\n"
	<< "  Source section:\n"
	<< "    space dim: " << field._spaceDim << "\n"
	<< "    vector field type: " << field._vecFieldType << "\n"
	<< "    scale: " << field._scale << "\n"
	<< "    size: " << srcSize
	<< "  Destination section:\n"
	<< "    space dim: " << _spaceDim << "\n"
	<< "    vector field type: " << _vecFieldType << "\n"
	<< "    scale: " << _scale << "\n"
	<< "    size: " << dstSize;
    throw std::runtime_error(msg.str());
  } // if

  // Add values from field
  const SieveRealSection::chart_type& chart = _section->getChart();
  const SieveRealSection::chart_type::const_iterator chartEnd = chart.end();

  // Assume fiber dimension is uniform
  const int fiberDim = _section->getFiberDimension(*chart.begin());
  double_array values(fiberDim);

  for (SieveRealSection::chart_type::const_iterator c_iter = chart.begin();
       c_iter != chartEnd;
       ++c_iter) {
    assert(fiberDim == field._section->getFiberDimension(*c_iter));
    assert(fiberDim == _section->getFiberDimension(*c_iter));
    field._section->restrictPoint(*c_iter, &values[0], values.size());
    _section->updateAddPoint(*c_iter, &values[0]);
  } // for
} // operator+=

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

  assert(!_section.isNull());
  const SieveRealSection::chart_type& chart = _section->getChart();
  const SieveRealSection::chart_type::const_iterator chartEnd = chart.end();

  // Assume fiber dimension is uniform
  const int fiberDim = _section->getFiberDimension(*chart.begin());
  double_array values(fiberDim);

  spatialdata::units::Nondimensional normalizer;

  for (SieveRealSection::chart_type::const_iterator c_iter = chart.begin();
       c_iter != chartEnd;
       ++c_iter) {
    assert(fiberDim == _section->getFiberDimension(*c_iter));

    _section->restrictPoint(*c_iter, &values[0], values.size());
    normalizer.dimensionalize(&values[0], values.size(), _scale);
    _section->updatePoint(*c_iter, &values[0]);
  } // for
  
} // dimensionalize


// End of file 
