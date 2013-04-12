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
// Copyright (c) 2010-2012 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

#include <portinfo>

#include "pylith/feassemble/Quadrature.hh" // USES Quadrature

#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/VisitorMesh.hh" // USES VecVisitorMesh

#include "pylith/utils/petscerror.h" // USES PYLITH_METHOD_BEGIN/END

// ----------------------------------------------------------------------
// Constructor
template<typename mesh_type, typename field_type>
pylith::meshio::CellFilterAvg<mesh_type, field_type>::CellFilterAvg(void) :
  _fieldAvg(0)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor
template<typename mesh_type, typename field_type>
pylith::meshio::CellFilterAvg<mesh_type, field_type>::~CellFilterAvg(void)
{ // destructor
  deallocate();
} // destructor  

// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
template<typename mesh_type, typename field_type>
void
pylith::meshio::CellFilterAvg<mesh_type, field_type>::deallocate(void)
{ // deallocate
  CellFilter<mesh_type, field_type>::deallocate();

  delete _fieldAvg; _fieldAvg = 0;
} // deallocate
  
// ----------------------------------------------------------------------
// Copy constructor.
template<typename mesh_type, typename field_type>
pylith::meshio::CellFilterAvg<mesh_type, field_type>::CellFilterAvg(
					       const CellFilterAvg& f) :
  CellFilter<mesh_type, field_type>(f),
  _fieldAvg(0)
{ // copy constructor
} // copy constructor

// ----------------------------------------------------------------------
// Create copy of filter.
template<typename mesh_type, typename field_type>
pylith::meshio::CellFilter<mesh_type, field_type>*
pylith::meshio::CellFilterAvg<mesh_type, field_type>::clone(void) const
{ // clone
  PYLITH_METHOD_BEGIN;
  
  pylith::meshio::CellFilter<mesh_type, field_type>* field = new CellFilterAvg<mesh_type,field_type>(*this);

  PYLITH_METHOD_RETURN(field);
} // clone

// ----------------------------------------------------------------------
// Get averaged field buffer.
template<typename mesh_type, typename field_type>
const field_type*
pylith::meshio::CellFilterAvg<mesh_type, field_type>::fieldAvg(void) const
{ // fieldAvg
  return _fieldAvg;
} // fieldAvg
  
// ----------------------------------------------------------------------
// Filter field.
template<typename mesh_type, typename field_type>
field_type&
pylith::meshio::CellFilterAvg<mesh_type,field_type>::filter(const field_type& fieldIn,
							    const char* label,
							    const int labelId)
{ // filter
  PYLITH_METHOD_BEGIN;

  const feassemble::Quadrature<mesh_type>* quadrature = CellFilter<mesh_type, field_type>::_quadrature;assert(quadrature);

  const int numQuadPts = quadrature->numQuadPts();
  const scalar_array& wts = quadrature->quadWts();
  
  PetscDM dmMesh = fieldIn.mesh().dmMesh();assert(dmMesh);
  PetscIS cellIS = NULL;
  PetscInt cStart, cEnd, numCells;
  PetscErrorCode err;

  if (!label) {
    PetscInt cMax;
    err = DMPlexGetHeightStratum(dmMesh, 0, &cStart, &cEnd);CHECK_PETSC_ERROR(err);
    err = DMPlexGetHybridBounds(dmMesh, &cMax, PETSC_NULL, PETSC_NULL, PETSC_NULL);CHECK_PETSC_ERROR(err);
    if (cMax >= 0) {cEnd = PetscMin(cEnd, cMax);}
    numCells = cEnd - cStart;
  } else {
    const PetscInt *cells;
    err = DMPlexGetStratumIS(dmMesh, label, 1, &cellIS);CHECK_PETSC_ERROR(err);
    err = ISGetSize(cellIS, &numCells);CHECK_PETSC_ERROR(err);
    err = ISGetIndices(cellIS, &cells);CHECK_PETSC_ERROR(err);
    cStart = cells[0];
    err = ISRestoreIndices(cellIS, &cells);CHECK_PETSC_ERROR(err);
    err = ISDestroy(&cellIS);CHECK_PETSC_ERROR(err);
  } // if

  topology::VecVisitorMesh fieldInVisitor(fieldIn);
  const PetscScalar* fieldInArray = fieldInVisitor.localArray();
  
  // Only processors with cells for output get the correct fiber dimension.
  PetscInt totalFiberDim = (numCells > 0) ? fieldInVisitor.sectionDof(cStart) : 0;
  const int fiberDim = totalFiberDim / numQuadPts;
  assert(fiberDim * numQuadPts == totalFiberDim);
  // Allocate field if necessary
  if (!_fieldAvg) {
    _fieldAvg = new field_type(fieldIn.mesh());assert(_fieldAvg);
    _fieldAvg->newSection(fieldIn, fiberDim);
    _fieldAvg->allocate();
  } else if (_fieldAvg->sectionSize() != numCells*fiberDim) {
    _fieldAvg->newSection(fieldIn, fiberDim);
    _fieldAvg->allocate();
  } // else

  assert(_fieldAvg);
  switch (fieldIn.vectorFieldType())
    { // switch
    case topology::FieldBase::MULTI_SCALAR:
      _fieldAvg->vectorFieldType(topology::FieldBase::SCALAR);
      break;
    case topology::FieldBase::MULTI_VECTOR:
      _fieldAvg->vectorFieldType(topology::FieldBase::VECTOR);
      break;
    case topology::FieldBase::MULTI_TENSOR:
      _fieldAvg->vectorFieldType(topology::FieldBase::TENSOR);
      break;
    case topology::FieldBase::MULTI_OTHER:
      _fieldAvg->vectorFieldType(topology::FieldBase::OTHER);
      break;
    case topology::FieldBase::SCALAR:
    case topology::FieldBase::VECTOR:
    case topology::FieldBase::TENSOR:
    case topology::FieldBase::OTHER:
    default :
      std::cerr << "Bad vector field type '" << fieldIn.vectorFieldType()
		<< "'." << std::endl;
      assert(0);
      throw std::logic_error("Bad vector field type for CellFilterAvg.");
    } // switch

  _fieldAvg->label(fieldIn.label());
  _fieldAvg->scale(fieldIn.scale());
  _fieldAvg->addDimensionOkay(true);

  topology::VecVisitorMesh fieldAvgVisitor(*_fieldAvg);
  PetscScalar* fieldAvgArray = fieldAvgVisitor.localArray();
  
  PylithScalar volume = 0.0;
  for (int iQuad=0; iQuad < numQuadPts; ++iQuad)
    volume += wts[iQuad];

  // Loop over cells
  if (cellIS) {
    const PetscInt *cells = NULL;
    err = ISGetIndices(cellIS, &cells);CHECK_PETSC_ERROR(err);
    for(PetscInt c = 0; c < numCells; ++c) {
      const PetscInt ioff = fieldInVisitor.sectionOffset(cells[c]);
      assert(totalFiberDim == fieldInVisitor.sectionDof(cells[c]));

      const PetscInt aoff = fieldAvgVisitor.sectionOffset(cells[c]);
      assert(fiberDim == fieldAvgVisitor.sectionDof(cells[c]));

      for(int i = 0; i < fiberDim; ++i) {
        fieldAvgArray[aoff+i] = 0.0;
        for(int iQuad = 0; iQuad < numQuadPts; ++iQuad)
          fieldAvgArray[aoff+i] += wts[iQuad] / volume * fieldInArray[ioff+iQuad*fiberDim+i];
      } // for
    } // for
    err = ISRestoreIndices(cellIS, &cells);CHECK_PETSC_ERROR(err);
    err = ISDestroy(&cellIS);CHECK_PETSC_ERROR(err);
  } else {
    for(PetscInt c = cStart; c < cEnd; ++c) {
      const PetscInt ioff = fieldInVisitor.sectionOffset(c);
      assert(totalFiberDim == fieldInVisitor.sectionDof(c));

      const PetscInt aoff = fieldAvgVisitor.sectionOffset(c);
      assert(fiberDim == fieldAvgVisitor.sectionDof(c));

      for(int i = 0; i < fiberDim; ++i) {
        fieldAvgArray[aoff+i] = 0.0;
        for(int iQuad = 0; iQuad < numQuadPts; ++iQuad)
          fieldAvgArray[aoff+i] += wts[iQuad] / volume * fieldInArray[ioff+iQuad*fiberDim+i];
      } // for
    } // for
  } // if/else
  PetscLogFlops(numCells * numQuadPts*fiberDim*3);

  PYLITH_METHOD_RETURN(*_fieldAvg);
} // filter


// End of file
