// -*- C++ -*-
//
// ======================================================================
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
// ======================================================================
//

#include <portinfo>

#include "ReverseCuthillMcKee.hh" // implementation of class methods

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/utils/error.h" // USES PYLITH_CHECK_ERROR

// ----------------------------------------------------------------------
// Reorder vertices and cells in mesh.
void
pylith::topology::ReverseCuthillMcKee::reorder(topology::Mesh* mesh)
{ // reorder
  assert(mesh);
  DMLabel label;
  PetscIS permutation;
  PetscDM dmOrig = mesh->dmMesh();
  PetscDM dmNew = NULL;
  PetscErrorCode err;

  err = DMPlexGetLabel(dmOrig, "material-id", &label);PYLITH_CHECK_ERROR(err);
  err = DMPlexGetOrdering(dmOrig, MATORDERINGRCM, label, &permutation);PYLITH_CHECK_ERROR(err);
  err = DMPlexPermute(dmOrig, permutation, &dmNew);PYLITH_CHECK_ERROR(err);
  err = ISDestroy(&permutation);PYLITH_CHECK_ERROR(err);
  mesh->dmMesh(dmNew);
} // reorder


// End of file 
