// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2011 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "BoundaryConditionPoints.hh" // implementation of object methods

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/FieldsNew.hh" // HOLDSA FieldsNew
#include "pylith/topology/Field.hh" // USES Field

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
  deallocate();
} // destructor

// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void 
pylith::bc::BoundaryConditionPoints::deallocate(void)
{ // deallocate
  BoundaryCondition::deallocate();
  delete _parameters; _parameters = 0;
} // deallocate
  
// ----------------------------------------------------------------------
// Get parameter fields.
const pylith::topology::FieldsNew<pylith::topology::Mesh>*
pylith::bc::BoundaryConditionPoints::parameterFields(void) const
{ // parameterFields
  return _parameters;
} // paramegetFields

// ----------------------------------------------------------------------
// Get mesh labels for points associated with boundary condition.
void
pylith::bc::BoundaryConditionPoints::_getPoints(const topology::Mesh& mesh)
{ // _getPoints
  typedef topology::Mesh::IntSection::chart_type chart_type;

  ALE::MemoryLogger& logger = ALE::MemoryLogger::singleton();
  logger.stagePush("BoundaryConditions");

  const ALE::Obj<SieveMesh>& sieveMesh = mesh.sieveMesh();
  assert(!sieveMesh.isNull());

  if (!sieveMesh->hasIntSection(_label)) {
    std::ostringstream msg;
    msg << "Could not find group of points '" << _label << "' in mesh.";
    throw std::runtime_error(msg.str());
  } // if

  const ALE::Obj<SieveMesh::int_section_type>& groupField = 
    sieveMesh->getIntSection(_label);
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

  logger.stagePop();
} // _getPoints


// End of file 
