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
// Copyright (c) 2010-2017 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

#include <portinfo>

#include "VertexFilterDecimateP1.hh" // Implementation of class methods

#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/Stratum.hh" // USES Stratum
#include "pylith/topology/VisitorMesh.hh" \
    // USES VecVisitorMesh

// ----------------------------------------------------------------------
// Constructor
pylith::meshio::VertexFilterDecimateP1::VertexFilterDecimateP1(void) :
    _fieldP1(0)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor
pylith::meshio::VertexFilterDecimateP1::~VertexFilterDecimateP1(void)
{ // destructor
    deallocate();
} // destructor

// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::meshio::VertexFilterDecimateP1::deallocate(void)
{ // deallocate
    PYLITH_METHOD_BEGIN;

    VertexFilter::deallocate();

    delete _fieldP1; _fieldP1 = 0;

    PYLITH_METHOD_END;
} // deallocate

// ----------------------------------------------------------------------
// Copy constructor.
pylith::meshio::VertexFilterDecimateP1::VertexFilterDecimateP1(const VertexFilterDecimateP1& f) :
    VertexFilter(f),
    _fieldP1(0)
{ // copy constructor
} // copy constructor

// ----------------------------------------------------------------------
// Create copy of filter.
pylith::meshio::VertexFilter*
pylith::meshio::VertexFilterDecimateP1::clone(void) const
{ // clone
    return new VertexFilterDecimateP1(*this);
} // clone

// ----------------------------------------------------------------------
// Filter field.
pylith::topology::Field&
pylith::meshio::VertexFilterDecimateP1::filter(const topology::Field& fieldIn)
{ // filter
    PYLITH_METHOD_BEGIN;

    PetscDM dmMesh = fieldIn.mesh().dmMesh();assert(dmMesh);
    pylith::topology::Stratum verticesStratum(dmMesh, pylith::topology::Stratum::DEPTH, 0);
    const PylithInt vStart = verticesStratum.begin();
    const PylithInt vEnd = verticesStratum.end();

    pylith::topology::VecVisitorMesh fieldInVisitor(fieldIn);

    // Allocate field if necessary
    if (!_fieldP1) {
        _fieldP1 = new pylith::topology::Field(fieldIn.mesh());
        _fieldP1->label(fieldIn.label());
        _fieldP1->dimensionalizeOkay(true);

        const pylith::string_vector& subfieldNames = fieldIn.subfieldNames();
        assert(1 == subfieldNames.size());
        const pylith::topology::Field::SubfieldInfo& info = fieldIn.subfieldInfo(subfieldNames[0].c_str());
        assert(0 == info.index);

        // Set discretization to P1
        pylith::topology::Field::Discretization feP1;
        feP1.basisOrder = 1;
        feP1.quadOrder = 1;
        feP1.isBasisContinuous = info.fe.isBasisContinuous;
        feP1.feSpace = info.fe.feSpace;

        _fieldP1->subfieldAdd(info.description, feP1);
        _fieldP1->subfieldsSetup();
        _fieldP1->allocate();
    } // if

    const PetscScalar* fieldInArray = fieldInVisitor.localArray();

    pylith::topology::VecVisitorMesh fieldP1Visitor(*_fieldP1);
    PetscScalar* fieldP1Array = fieldP1Visitor.localArray();

    // Loop over vertices
    for (PylithInt v = vStart; v < vEnd; ++v) {
        const PylithInt offIn = fieldInVisitor.sectionOffset(v);
        const PylithInt ndof = fieldInVisitor.sectionDof(v);
        assert(fieldP1Visitor.sectionDof(v) == ndof);

        const PylithInt offP1 = fieldP1Visitor.sectionOffset(v);

        for (PylithInt d = 0; d < ndof; ++d) {
            fieldP1Array[offP1+d] = fieldInArray[offIn+d];
        } // for
    } // for

    PYLITH_METHOD_RETURN(*_fieldP1);
} // filter


// End of file
