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

#include "CellFilterAvg.hh" // Implementation of class methods

#include "pylith/feassemble/Quadrature.hh" // USES Quadrature

#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/VisitorMesh.hh" // USES VecVisitorMesh
#include "pylith/topology/Stratum.hh" // USES StratumIS

#include "pylith/utils/error.h" // USES PYLITH_METHOD_BEGIN/END

// ----------------------------------------------------------------------
// Constructor
pylith::meshio::CellFilterAvg::CellFilterAvg(void) :
  _fieldAvg(0)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor
pylith::meshio::CellFilterAvg::~CellFilterAvg(void)
{ // destructor
  deallocate();
} // destructor  

// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::meshio::CellFilterAvg::deallocate(void)
{ // deallocate
  CellFilter::deallocate();

  delete _fieldAvg; _fieldAvg = 0;
} // deallocate
  
// ----------------------------------------------------------------------
// Copy constructor.
pylith::meshio::CellFilterAvg::CellFilterAvg(const CellFilterAvg& f) :
  CellFilter(f),
  _fieldAvg(0)
{ // copy constructor
} // copy constructor

// ----------------------------------------------------------------------
// Create copy of filter.
pylith::meshio::CellFilter*
pylith::meshio::CellFilterAvg::clone(void) const
{ // clone
  PYLITH_METHOD_BEGIN;
  
  pylith::meshio::CellFilter* f = new CellFilterAvg(*this);

  PYLITH_METHOD_RETURN(f);
} // clone

// ----------------------------------------------------------------------
// Get averaged field buffer.
const pylith::topology::Field*
pylith::meshio::CellFilterAvg::fieldAvg(void) const
{ // fieldAvg
  return _fieldAvg;
} // fieldAvg
  
// ----------------------------------------------------------------------
// Filter field.
pylith::topology::Field&
pylith::meshio::CellFilterAvg::filter(const topology::Field& fieldIn,
				      const char* label,
				      const int labelId)
{ // filter
  PYLITH_METHOD_BEGIN;

  const feassemble::Quadrature* quadrature = CellFilter::_quadrature;assert(quadrature);

  const int numQuadPts = quadrature->numQuadPts();
  const scalar_array& wts = quadrature->quadWts();
  
  PetscDM dmMesh = fieldIn.mesh().dmMesh();assert(dmMesh);
  PetscInt cStart = 0, cEnd, numCells;
  PetscErrorCode err;

  if (!label) {
    PetscInt h, cMax;
    err = DMPlexGetVTKCellHeight(dmMesh, &h);PYLITH_CHECK_ERROR(err);
    err = DMPlexGetHeightStratum(dmMesh, h, &cStart, &cEnd);PYLITH_CHECK_ERROR(err);
    err = DMPlexGetHybridBounds(dmMesh, &cMax, NULL, NULL, NULL);PYLITH_CHECK_ERROR(err);
    if (cMax >= 0) {cEnd = PetscMin(cEnd, cMax);}
    numCells = cEnd - cStart;
  } else {
    if (!_cellsIS) {
      const bool includeOnlyCells = true;
      _cellsIS = new topology::StratumIS(dmMesh, label, labelId, includeOnlyCells);assert(_cellsIS);
    } // if
    numCells = _cellsIS->size();
    if (numCells > 0) {
      cStart = _cellsIS->points()[0];
    } // if
  } // if/else

  topology::VecVisitorMesh fieldInVisitor(fieldIn);
  const PetscScalar* fieldInArray = fieldInVisitor.localArray();
  
  // Only processors with cells for output get the correct fiber dimension.
  PetscInt totalFiberDim = (numCells > 0) ? fieldInVisitor.sectionDof(cStart) : 0;
  const int fiberDim = totalFiberDim / numQuadPts;
  assert(fiberDim * numQuadPts == totalFiberDim);
  // The decision to reallocate a field must be collective
  PetscInt reallocate = ((!_fieldAvg) || (_fieldAvg->sectionSize() != numCells*fiberDim));
  PetscInt reallocateGlobal = 0;
  err = MPI_Allreduce(&reallocate, &reallocateGlobal, 1, MPIU_INT, MPI_LOR, fieldIn.mesh().comm());PYLITH_CHECK_ERROR(err);
  if (reallocateGlobal) {
    if (!_fieldAvg) {
      _fieldAvg = new topology::Field(fieldIn.mesh());assert(_fieldAvg);
    } // if
    _fieldAvg->newSection(fieldIn, fiberDim);
    _fieldAvg->allocate();
  } // if

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
      std::ostringstream msg;
      msg << "Bad vector field type '" << fieldIn.vectorFieldType() << " for CellFilterAvg'." << std::endl;
      throw std::logic_error(msg.str());
    } // switch

  _fieldAvg->label(fieldIn.label());
  _fieldAvg->scale(fieldIn.scale());
  _fieldAvg->dimensionalizeOkay(true);

  topology::VecVisitorMesh fieldAvgVisitor(*_fieldAvg);
  PetscScalar* fieldAvgArray = fieldAvgVisitor.localArray();
  
  PylithScalar volume = 0.0;
  for (int iQuad=0; iQuad < numQuadPts; ++iQuad)
    volume += wts[iQuad];

  // Loop over cells
  if (_cellsIS) {
    const PetscInt* cells = _cellsIS->points();
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
