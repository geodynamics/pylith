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
// Copyright (c) 2010-2016 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "DirichletNew.hh" // implementation of object methods

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/VisitorMesh.hh" // USES VecVisitorMesh
#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include <strings.h> // USES strcasecmp()
#include <cassert> // USES assert()
#include <stdexcept> // USES std::runtime_error
#include <sstream> // USES std::ostringstream

// ----------------------------------------------------------------------
// Default constructor.
pylith::bc::DirichletNew::DirichletNew(void) :
    _boundaryMesh(0)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::bc::DirichletNew::~DirichletNew(void)
{ // destructor
    deallocate();
} // destructor

// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::bc::DirichletNew::deallocate(void)
{ // deallocate
    PYLITH_METHOD_BEGIN;

    ConstraintPointwise::deallocate();

    delete _boundaryMesh; _boundaryMesh = 0;

    PYLITH_METHOD_END;
} // deallocate

// ----------------------------------------------------------------------
// Initialize boundary condition.
void
pylith::bc::DirichletNew::initialize(const pylith::topology::Field& solution)
{ // initialize
    PYLITH_METHOD_BEGIN;

    _boundaryMesh = new topology::Mesh(solution.mesh(), _label.c_str()); assert(_boundaryMesh);

    PYLITH_METHOD_END;
} // initialize


// End of file
