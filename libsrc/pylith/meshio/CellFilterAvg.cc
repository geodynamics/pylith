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
  return new CellFilterAvg<mesh_type,field_type>(*this);
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
pylith::meshio::CellFilterAvg<mesh_type,field_type>::filter(
						const field_type& fieldIn,
						const char* label,
						const int labelId)
{ // filter
  typedef typename mesh_type::SieveMesh SieveMesh;
  typedef typename SieveMesh::label_sequence label_sequence;
  typedef typename field_type::Mesh::RealSection RealSection;

  const feassemble::Quadrature<mesh_type>* quadrature = 
    CellFilter<mesh_type, field_type>::_quadrature;
  assert(0 != quadrature);

  const int numQuadPts = quadrature->numQuadPts();
  const scalar_array& wts = quadrature->quadWts();
  
  DM             dmMesh = fieldIn.mesh().dmMesh();
  IS             cellIS = PETSC_NULL;
  PetscInt       cStart, cEnd, numCells;
  PetscErrorCode err;

  assert(dmMesh);
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
  }

  // Only processors with cells for output get the correct fiber dimension.
  PetscSection sectionIn     = fieldIn.petscSection();
  Vec          vecIn         = fieldIn.localVector();
  PetscInt     totalFiberDim = 0;

  assert(sectionIn);
  if (numCells) {err = PetscSectionGetDof(sectionIn, cStart, &totalFiberDim);CHECK_PETSC_ERROR(err);}
  const int fiberDim = totalFiberDim / numQuadPts;
  assert(fiberDim * numQuadPts == totalFiberDim);
  // Allocate field if necessary
  ALE::MemoryLogger& logger = ALE::MemoryLogger::singleton();
  logger.stagePush("OutputFields");
  if (0 == _fieldAvg) {
    _fieldAvg = new field_type(fieldIn.mesh());
    assert(0 != _fieldAvg);
    _fieldAvg->newSection(fieldIn, fiberDim);
    _fieldAvg->allocate();
  } else if (_fieldAvg->sectionSize() != numCells*fiberDim) {
    _fieldAvg->newSection(fieldIn, fiberDim);
    _fieldAvg->allocate();
  } // else
  logger.stagePop();

  assert(0 != _fieldAvg);
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
  PetscSection sectionAvg = _fieldAvg->petscSection();
  Vec          vecAvg     = _fieldAvg->localVector();
  _fieldAvg->label(fieldIn.label());
  _fieldAvg->scale(fieldIn.scale());
  _fieldAvg->addDimensionOkay(true);

  PylithScalar volume = 0.0;
  for (int iQuad=0; iQuad < numQuadPts; ++iQuad)
    volume += wts[iQuad];

  // Loop over cells
  PetscScalar *arrayIn, *arrayAvg;

  err = VecGetArray(vecIn,  &arrayIn);CHECK_PETSC_ERROR(err);
  err = VecGetArray(vecAvg, &arrayAvg);CHECK_PETSC_ERROR(err);
  if (cellIS) {
    const PetscInt *cells;

    err = ISGetIndices(cellIS, &cells);CHECK_PETSC_ERROR(err);
    for(PetscInt c = 0; c < numCells; ++c) {
      PetscInt dof, off, adof, aoff;
    
      err = PetscSectionGetDof(sectionIn, cells[c], &dof);CHECK_PETSC_ERROR(err);
      err = PetscSectionGetOffset(sectionIn, cells[c], &off);CHECK_PETSC_ERROR(err);
      err = PetscSectionGetDof(sectionAvg, cells[c], &adof);CHECK_PETSC_ERROR(err);
      err = PetscSectionGetOffset(sectionAvg, cells[c], &aoff);CHECK_PETSC_ERROR(err);
      assert(totalFiberDim == dof);
      assert(fiberDim == adof);
      for(int i = 0; i < adof; ++i) {
        arrayAvg[aoff+i] = 0.0;
        for(int iQuad = 0; iQuad < numQuadPts; ++iQuad)
          arrayAvg[aoff+i] += wts[iQuad] / volume * arrayIn[off+iQuad*fiberDim+i];
      }
    }
    err = ISRestoreIndices(cellIS, &cells);CHECK_PETSC_ERROR(err);
    err = ISDestroy(&cellIS);CHECK_PETSC_ERROR(err);
  } else {
    for(PetscInt c = cStart; c < cEnd; ++c) {
      PetscInt dof, off, adof, aoff;
    
      err = PetscSectionGetDof(sectionIn, c, &dof);CHECK_PETSC_ERROR(err);
      err = PetscSectionGetOffset(sectionIn, c, &off);CHECK_PETSC_ERROR(err);
      err = PetscSectionGetDof(sectionAvg, c, &adof);CHECK_PETSC_ERROR(err);
      err = PetscSectionGetOffset(sectionAvg, c, &aoff);CHECK_PETSC_ERROR(err);
      assert(totalFiberDim == dof);
      assert(fiberDim == adof);
      for(int i = 0; i < adof; ++i) {
        arrayAvg[aoff+i] = 0.0;
        for(int iQuad = 0; iQuad < numQuadPts; ++iQuad)
          arrayAvg[aoff+i] += wts[iQuad] / volume * arrayIn[off+iQuad*fiberDim+i];
      }
    } // for
  }
  err = VecRestoreArray(vecIn, &arrayIn);CHECK_PETSC_ERROR(err);
  err = VecRestoreArray(vecAvg, &arrayAvg);CHECK_PETSC_ERROR(err);
  PetscLogFlops(numCells * numQuadPts*fiberDim*3);

  return *_fieldAvg;
} // filter


// End of file
