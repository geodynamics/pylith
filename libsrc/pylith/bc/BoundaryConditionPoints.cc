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
// Copyright (c) 2010-2015 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "BoundaryConditionPoints.hh" // implementation of object methods

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Fields.hh" // HOLDSA Fields
#include "pylith/topology/Field.hh" // USES Field

#include <stdexcept> // USES std::runtime_error()

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
  PYLITH_METHOD_BEGIN;

  BoundaryCondition::deallocate();
  delete _parameters; _parameters = 0;

  PYLITH_METHOD_END;
} // deallocate
  
// ----------------------------------------------------------------------
// Get parameter fields.
const pylith::topology::Fields*
pylith::bc::BoundaryConditionPoints::parameterFields(void) const
{ // parameterFields
  return _parameters;
} // paramegetFields

// ----------------------------------------------------------------------
// Get mesh labels for points associated with boundary condition.
void
pylith::bc::BoundaryConditionPoints::_getPoints(const topology::Mesh& mesh)
{ // _getPoints
  PYLITH_METHOD_BEGIN;

  PetscDM dmMesh = mesh.dmMesh();assert(dmMesh);
  PetscDMLabel label = NULL;
  PetscIS pointIS = NULL;
  const PetscInt *points;
  PetscInt numPoints = 0, vStart = 0, vEnd = 0, numVertices = 0;
  PetscBool hasLabel = PETSC_FALSE;
  PetscErrorCode err;
  err = DMHasLabel(dmMesh, _label.c_str(), &hasLabel);PYLITH_CHECK_ERROR(err);
  if (!hasLabel) {
    std::ostringstream msg;
    msg << "Could not find group of points '" << _label << "' in mesh.";
    throw std::runtime_error(msg.str());
  } // if

  err = DMPlexGetDepthStratum(dmMesh, 0, &vStart, &vEnd);PYLITH_CHECK_ERROR(err);
  err = DMGetLabel(dmMesh, _label.c_str(), &label);PYLITH_CHECK_ERROR(err);
  err = DMLabelGetStratumIS(label, 1, &pointIS);PYLITH_CHECK_ERROR(err);
  if (pointIS) {
    err = ISGetLocalSize(pointIS, &numPoints);PYLITH_CHECK_ERROR(err);
    err = ISGetIndices(pointIS, &points);PYLITH_CHECK_ERROR(err);
  } // if

  // For now, only use vertices
  for (PetscInt p = 0; p < numPoints; ++p) {
    if ((points[p] >= vStart) && (points[p] < vEnd)) {
      ++numVertices;
    } // if
  } // for
  _points.resize(numVertices);
  for (PetscInt p = 0, v = 0; p < numPoints; ++p) {
    if ((points[p] >= vStart) && (points[p] < vEnd)) {
      _points[v++] = points[p];
    } // if
  } // for
  if (pointIS) {
    err = ISRestoreIndices(pointIS, &points);PYLITH_CHECK_ERROR(err);
    err = ISDestroy(&pointIS);PYLITH_CHECK_ERROR(err);
  } // if

  PYLITH_METHOD_END;
} // _getPoints


// End of file 
