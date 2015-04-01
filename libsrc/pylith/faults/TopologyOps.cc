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

#include "TopologyOps.hh" // implementation of object methods

#include <cassert> // USES assert()

// ----------------------------------------------------------------------
void
pylith::faults::TopologyOps::classifyCellsDM(PetscDM dmMesh,
                                             PetscInt vertex,
                                             const int depth,
                                             const int faceSize,
                                             PetscInt firstCohesiveCell,
                                             PointSet& replaceCells,
                                             PointSet& noReplaceCells,
                                             const int debug)
{
  // Replace all cells on a given side of the fault with a vertex on the fault
  PointSet        vReplaceCells;
  PointSet        vNoReplaceCells;
  const PetscInt *support;
  PetscInt        supportSize, s, classifyTotal = 0;
  PetscBool       modified = PETSC_FALSE;
  PetscErrorCode  err;

  if (debug) {std::cout << "Checking fault vertex " << vertex << std::endl;}
  err = DMPlexGetSupportSize(dmMesh, vertex, &supportSize);PYLITH_CHECK_ERROR(err);
  err = DMPlexGetSupport(dmMesh, vertex, &support);PYLITH_CHECK_ERROR(err);
  for (s = 0; s < supportSize; ++s) {
    const PetscInt point = support[s];

    if (point >= firstCohesiveCell) return;
    if (replaceCells.find(point)   != replaceCells.end())   vReplaceCells.insert(point);
    if (noReplaceCells.find(point) != noReplaceCells.end()) vNoReplaceCells.insert(point);
    modified = PETSC_TRUE;
    ++classifyTotal;
  }
  PetscInt classifySize = vReplaceCells.size() + vNoReplaceCells.size();

  while (modified && (classifySize < classifyTotal)) {
    modified = PETSC_FALSE;
    for (s = 0; s < supportSize; ++s) {
      const PetscInt point      = support[s];
      PetscBool      classified = PETSC_FALSE;
    
      if (debug) {
        const PetscInt *cone;
        PetscInt        coneSize;

        std::cout << "Checking neighbor " << vertex << std::endl;
        err = DMPlexGetConeSize(dmMesh, vertex, &coneSize);PYLITH_CHECK_ERROR(err);
        err = DMPlexGetCone(dmMesh, vertex, &cone);PYLITH_CHECK_ERROR(err);
        for (PetscInt c = 0; c < coneSize; ++c) {
          std::cout << "  cone point " << cone[c] << std::endl;
        }
      }
      if (vReplaceCells.find(point) != vReplaceCells.end()) {
        if (debug) std::cout << "  already in replaceCells" << std::endl;
        continue;
      } // if
      if (vNoReplaceCells.find(point) != vNoReplaceCells.end()) {
        if (debug) std::cout << "  already in noReplaceCells" << std::endl;
        continue;
      } // if
      if (point >= firstCohesiveCell) {
        if (debug) std::cout << "  already a cohesive cell" << std::endl;
        continue;
      } // if
      // If neighbor shares a face with anyone in replaceCells, then add
      for (PointSet::const_iterator c_iter = vReplaceCells.begin(); c_iter != vReplaceCells.end(); ++c_iter) {
        const PetscInt *coveringPoints;
        PetscInt        numCoveringPoints, points[2];

        points[0] = point; points[1] = *c_iter;
        err = DMPlexGetMeet(dmMesh, 2, points, &numCoveringPoints, &coveringPoints);PYLITH_CHECK_ERROR(err);
        err = DMPlexRestoreMeet(dmMesh, 2, points, &numCoveringPoints, &coveringPoints);PYLITH_CHECK_ERROR(err);
        if (numCoveringPoints == faceSize) {
          if (debug) std::cout << "    Scheduling " << point << " for replacement" << std::endl;
          vReplaceCells.insert(point);
          modified   = PETSC_TRUE;
          classified = PETSC_TRUE;
          break;
        } // if
      } // for
      if (classified) continue;
      // It is unclear whether taking out the noReplace cells will speed this up
      for (PointSet::const_iterator c_iter = vNoReplaceCells.begin(); c_iter != vNoReplaceCells.end(); ++c_iter) {
        const PetscInt *coveringPoints;
        PetscInt        numCoveringPoints, points[2];

        points[0] = point; points[1] = *c_iter;
        err = DMPlexGetMeet(dmMesh, 2, points, &numCoveringPoints, &coveringPoints);PYLITH_CHECK_ERROR(err);
        err = DMPlexRestoreMeet(dmMesh, 2, points, &numCoveringPoints, &coveringPoints);PYLITH_CHECK_ERROR(err);
        if (numCoveringPoints == faceSize) {
          if (debug) std::cout << "    Scheduling " << point << " for no replacement" << std::endl;
          vNoReplaceCells.insert(point);
          modified   = PETSC_TRUE;
          classified = PETSC_TRUE;
          break;
        } // for
      } // for
    }
    if (debug) {
      std::cout << "classifySize: " << classifySize << std::endl;
      std::cout << "classifyTotal: " << classifyTotal << std::endl;
      std::cout << "vReplaceCells.size: " << vReplaceCells.size() << std::endl;
      std::cout << "vNoReplaceCells.size: " << vNoReplaceCells.size() << std::endl;
    }
    assert(size_t(classifySize) < vReplaceCells.size() + vNoReplaceCells.size());
    classifySize = vReplaceCells.size() + vNoReplaceCells.size();
    if (classifySize > classifyTotal) {
      std::ostringstream msg;
      msg << "Internal error classifying cells during creation of cohesive cells."
          << "  classifySize: " << classifySize << ", classifyTotal: " << classifyTotal;
      throw std::logic_error(msg.str());
    } // if
  }
  replaceCells.insert(vReplaceCells.begin(), vReplaceCells.end());
  // More checking
  noReplaceCells.insert(vNoReplaceCells.begin(), vNoReplaceCells.end());
}

// End of file
