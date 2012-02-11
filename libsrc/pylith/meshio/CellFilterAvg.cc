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
// Copyright (c) 2010-2011 University of California, Davis
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
  
  const ALE::Obj<SieveMesh>& sieveMesh = fieldIn.mesh().sieveMesh();
  assert(!sieveMesh.isNull());
  const int cellDepth = (sieveMesh->depth() == -1) ? -1 : 1;
  const int depth = (0 == label) ? cellDepth : labelId;
  const std::string labelName = (0 == label) ?
    ((sieveMesh->hasLabel("censored depth")) ?
     "censored depth" : "depth") : label;

  const ALE::Obj<label_sequence>& cells = 
    sieveMesh->getLabelStratum(labelName, depth);
  assert(!cells.isNull());
  const typename label_sequence::iterator cellsBegin = cells->begin();
  const typename label_sequence::iterator cellsEnd = cells->end();

  // Only processors with cells for output get the correct fiber dimension.
  const ALE::Obj<RealSection>& sectionIn = fieldIn.section();
  assert(!sectionIn.isNull());
  const int totalFiberDim = (cellsBegin != cellsEnd) ?
    sectionIn->getFiberDimension(*cellsBegin) : 0;
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
  } else if (_fieldAvg->sectionSize() != cells->size()*fiberDim) {
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
  const ALE::Obj<RealSection>& sectionAvg = _fieldAvg->section();
  _fieldAvg->label(fieldIn.label());
  _fieldAvg->scale(fieldIn.scale());
  _fieldAvg->addDimensionOkay(true);

  scalar_array fieldAvgCell(fiberDim);
  PylithScalar scalar = 0.0;
  for (int iQuad=0; iQuad < numQuadPts; ++iQuad)
    scalar += wts[iQuad];

  // Loop over cells
  for (typename label_sequence::iterator c_iter=cellsBegin;
       c_iter != cellsEnd;
       ++c_iter) {
    const PylithScalar* values = sectionIn->restrictPoint(*c_iter);
    assert(totalFiberDim == sectionIn->getFiberDimension(*c_iter));
    
    fieldAvgCell = 0.0;
    for (int iQuad=0; iQuad < numQuadPts; ++iQuad)
      for (int i=0; i < fiberDim; ++i)
	fieldAvgCell[i] += wts[iQuad] / scalar * values[iQuad*fiberDim+i];

    sectionAvg->updatePoint(*c_iter, &fieldAvgCell[0]);
  } // for
  PetscLogFlops( cells->size() * numQuadPts*fiberDim*3 );

  return *_fieldAvg;
} // filter


// End of file
