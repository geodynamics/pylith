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
  assert(!_section.isNull());

  const SieveMesh::point_type pointMin = 
    *std::min_element(points->begin(), points->end());
  const SieveMesh::point_type pointMax = 
    *std::max_element(points->begin(), points->end());
  _section->setChart(SieveRealSection::chart_type(pointMin, pointMax+1));
  _section->setFiberDimension(points, _fiberDim);
  _mesh->allocate(_section);
} // createSection

// ----------------------------------------------------------------------
// Create section given chart.
void
pylith::topology::FieldUniform::createSection(
			const SieveRealSection::chart_type& chart)
{ // createSection
  assert(!_section.isNull());

  _section->setChart(chart);

  const SieveRealSection::chart_type::const_iterator chartEnd = chart.end();
  for (SieveRealSection::chart_type::const_iterator c_iter = chart.begin();
       c_iter != chartEnd;
       ++c_iter)
    _section->setFiberDimension(*c_iter, _fiberDim);
  _mesh->allocate(_section);
} // createSection


// End of file 
