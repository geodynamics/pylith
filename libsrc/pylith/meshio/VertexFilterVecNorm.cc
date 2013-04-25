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
// Copyright (c) 2010-2013 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

#include <portinfo>

#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/Stratum.hh" // USES Stratum
#include "pylith/topology/VisitorMesh.hh" // USES VecVisitorMesh

// ----------------------------------------------------------------------
// Constructor
template<typename field_type>
pylith::meshio::VertexFilterVecNorm<field_type>::VertexFilterVecNorm(void) :
  _fieldVecNorm(0)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor
template<typename field_type>
pylith::meshio::VertexFilterVecNorm<field_type>::~VertexFilterVecNorm(void)
{ // destructor
  deallocate();
} // destructor  

// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
template<typename field_type>
void
pylith::meshio::VertexFilterVecNorm<field_type>::deallocate(void)
{ // deallocate
  PYLITH_METHOD_BEGIN;

  VertexFilter<field_type>::deallocate();  

  delete _fieldVecNorm; _fieldVecNorm = 0;

  PYLITH_METHOD_END;
} // deallocate
  
// ----------------------------------------------------------------------
// Copy constructor.
template<typename field_type>
pylith::meshio::VertexFilterVecNorm<field_type>::VertexFilterVecNorm(const VertexFilterVecNorm& f) :
  VertexFilter<field_type>(f),
  _fieldVecNorm(0)
{ // copy constructor
} // copy constructor

// ----------------------------------------------------------------------
// Create copy of filter.
template<typename field_type>
pylith::meshio::VertexFilter<field_type>*
pylith::meshio::VertexFilterVecNorm<field_type>::clone(void) const
{ // clone
  return new VertexFilterVecNorm(*this);
} // clone

// ----------------------------------------------------------------------
// Filter field.
template<typename field_type>
field_type&
pylith::meshio::VertexFilterVecNorm<field_type>::filter(const field_type& fieldIn)
{ // filter
  PYLITH_METHOD_BEGIN;

  PetscDM dmMesh = fieldIn.mesh().dmMesh();assert(dmMesh);
  topology::Stratum verticesStratum(dmMesh, topology::Stratum::DEPTH, 0);
  const PetscInt vStart = verticesStratum.begin();
  const PetscInt vEnd = verticesStratum.end();

  topology::VecVisitorMesh fieldInVisitor(fieldIn);

  // Only processors with cells for output get the correct fiber dimension.
  PetscInt fiberDimIn = (verticesStratum.size() > 0) ? fieldInVisitor.sectionDof(vStart) : 0;
  const int fiberDimNorm = 1;

  // Allocate field if necessary
  if (!_fieldVecNorm) {
    _fieldVecNorm = new field_type(fieldIn.mesh());
    _fieldVecNorm->label("vector norm");
    _fieldVecNorm->newSection(fieldIn, fiberDimNorm);
    _fieldVecNorm->allocate();

    _fieldVecNorm->label(fieldIn.label());
    _fieldVecNorm->scale(fieldIn.scale());
    switch (fieldIn.vectorFieldType())
      { // switch
      case topology::FieldBase::SCALAR:
      case topology::FieldBase::VECTOR:
	_fieldVecNorm->vectorFieldType(topology::FieldBase::SCALAR);
	break;
      case topology::FieldBase::MULTI_SCALAR:
      case topology::FieldBase::MULTI_VECTOR:
      case topology::FieldBase::MULTI_TENSOR:
      case topology::FieldBase::MULTI_OTHER:
      case topology::FieldBase::TENSOR:
      case topology::FieldBase::OTHER:
      default :
	std::cerr << "Bad vector field type '" << fieldIn.vectorFieldType()
		  << "'." << std::endl;
	assert(0);
	throw std::logic_error("Bad vector field type in VertexFilterVecNorm.");
      } // switch
  } // if

  const PetscScalar* fieldInArray = fieldInVisitor.localArray();

  topology::VecVisitorMesh fieldNormVisitor(*_fieldVecNorm);
  PetscScalar* fieldNormArray = fieldNormVisitor.localArray();

  // Loop over vertices
  for(PetscInt v = vStart; v < vEnd; ++v) {
    const PetscInt ioff = fieldInVisitor.sectionOffset(v);
    assert(fiberDimIn == fieldInVisitor.sectionDof(v));

    const PetscInt noff = fieldNormVisitor.sectionOffset(v);

    PylithScalar norm = 0.0;
    for(PetscInt d = 0; d < fiberDimIn; ++d) {
      norm += fieldInArray[ioff+d]*fieldInArray[ioff+d];
    } // for
    fieldNormArray[noff] = sqrt(norm);
  } // for
  PetscLogFlops((vEnd-vStart) * (1 + 2*fiberDimIn));

  PYLITH_METHOD_RETURN(*_fieldVecNorm);
} // filter


// End of file
