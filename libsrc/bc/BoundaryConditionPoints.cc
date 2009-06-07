// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
//                           Brad T. Aagaard
//                        U.S. Geological Survey
//
// {LicenseText}
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "BoundaryConditionPoints.hh" // implementation of object methods

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/Fields.hh" // USES Fields

#include <stdexcept> // USES std::runtime_error()

// ----------------------------------------------------------------------
typedef pylith::topology::Mesh::SieveMesh SieveMesh;

// ----------------------------------------------------------------------
// Default constructor.
pylith::bc::BoundaryConditionPoints::BoundaryConditionPoints(void) :
  _parameters(0)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::bc::BoundaryConditionPoints::~BoundaryConditionPoints(void)
{ // destructor
  delete _parameters; _parameters = 0;
} // destructor

// ----------------------------------------------------------------------
// Get mesh labels for points associated with boundary condition.
void
pylith::bc::BoundaryConditionPoints::_getPoints(const topology::Mesh& mesh)
{ // _getPoints
  typedef topology::Mesh::IntSection::chart_type chart_type;

  const ALE::Obj<SieveMesh>& sieveMesh = mesh.sieveMesh();
  assert(!sieveMesh.isNull());

  const ALE::Obj<SieveMesh::int_section_type>& groupField = 
    sieveMesh->getIntSection(_label);
  if (groupField.isNull()) {
    std::ostringstream msg;
    msg << "Could not find group of points '" << _label << "' in mesh.";
    throw std::runtime_error(msg.str());
  } // if
  assert(!groupField.isNull());
  const chart_type& chart = groupField->getChart();
  const chart_type::const_iterator& chartEnd = chart.end();
  const int numPoints = groupField->size();
  _points.resize(numPoints);
  int i = 0;
  for(chart_type::const_iterator c_iter = chart.begin();
      c_iter != chartEnd;
      ++c_iter)
    if (groupField->getFiberDimension(*c_iter))
      _points[i++] = *c_iter;
} // _getPoints


// End of file 
