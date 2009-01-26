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

#include "FieldUniform.hh" // implementation of class methods

#include <stdexcept> // USES std::runtime_error
#include <sstream> // USES std::ostringstream
#include <cassert> // USES assert()

// ----------------------------------------------------------------------
// Default constructor.
pylith::topology::FieldUniform::FieldUniform(const ALE::Obj<SieveMesh>& mesh,
					     const int fiberDim) :
  Field(mesh),
  _fiberDim(fiberDim)
{ // constructor
  assert(fiberDim >= 0);
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::topology::FieldUniform::~FieldUniform(void)
{ // destructor
} // destructor

// ----------------------------------------------------------------------
// Create section given points.
void
pylith::topology::FieldUniform::createSection(
			const ALE::Obj<SieveMesh::label_sequence>& points)
{ // createSection
  newSection(points, _fiberDim);
  _mesh->allocate(_section);
} // createSection

// ----------------------------------------------------------------------
// Create section given chart.
void
pylith::topology::FieldUniform::createSection(
			const SieveRealSection::chart_type& chart)
{ // createSection
  if (_section.isNull())
    newSection();

  _section->setChart(chart);

  const SieveRealSection::chart_type::const_iterator chartEnd = chart.end();
  for (SieveRealSection::chart_type::const_iterator c_iter = chart.begin();
       c_iter != chartEnd;
       ++c_iter)
    _section->setFiberDimension(*c_iter, _fiberDim);
  _mesh->allocate(_section);
} // createSection


// End of file 
