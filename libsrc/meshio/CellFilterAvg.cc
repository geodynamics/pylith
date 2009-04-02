// -*- C++ -*-
//
// ======================================================================
//
//                           Brad T. Aagaard
//                        U.S. Geological Survey
//
// {LicenseText}
//
// ======================================================================
//

#include <portinfo>

#include "pylith/feassemble/Quadrature.hh" // USES Quadrature

// ----------------------------------------------------------------------
// Constructor
template<typename mesh_type>
pylith::meshio::CellFilterAvg<mesh_type>::CellFilterAvg(void)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor
template<typename mesh_type>
pylith::meshio::CellFilterAvg<mesh_type>::~CellFilterAvg(void)
{ // destructor
} // destructor  

// ----------------------------------------------------------------------
// Copy constructor.
template<typename mesh_type>
pylith::meshio::CellFilterAvg<mesh_type>::CellFilterAvg(const CellFilterAvg& f) :
  CellFilter(f)
{ // copy constructor
} // copy constructor

// ----------------------------------------------------------------------
// Create copy of filter.
template<typename mesh_type>
pylith::meshio::CellFilter<mesh_type>*
pylith::meshio::CellFilterAvg<mesh_type>::clone(void) const
{ // clone
  return new CellFilterAvg<mesh_type>(*this);
} // clone

// ----------------------------------------------------------------------
// Filter field.
template<typename mesh_type>
const pylith::topology::Field<mesh_type>&
pylith::meshio::CellFilterAvg<mesh_type>::filter(
				  const topology::Field<mesh_type>& fieldIn,
				  const char* label,
				  const int labelId)
{ // filter
  assert(0 != _quadrature);

  const int numQuadPts = _quadrature->numQuadPts();
  const double_array& wts = _quadrature->quadWts();
  
  const ALE::Obj<SieveMesh>& sieveMesh = fieldIn.mesh().sieveMesh();
  assert(!sieveMesh.isNull());

  const ALE::Obj<SieveMesh::label_sequence>& cells = (0 == label) ?
    mesh->heightStratum(0) :
    mesh->getLabelStratum(label, labelId);
  assert(!cells.isNull());
  const SieveMesh::label_sequence::iterator cellsEnd = cells->end();

  // Only processors with cells for output get the correct fiber dimension.
  const ALE::Obj<RealSection>& sectionIn = fieldIn.section();
  assert(!sectionIn.isNull());
  const int totalFiberDim = sectionIn->getFiberDimension(*cells->begin());
  const int fiberDim = totalFiberDim / numQuadPts;
  assert(fiberDim * numQuadPts == totalFiberDim);

  // Allocate field if necessary
  if (0 == _fieldAvg) {
    _fieldAvg = new topology::Field<mesh_type>(fieldIn.mesh());
    assert(0 != _fieldAvg);
    _fieldAvg->newSection(fieldIn->getChart(), fiberDim);
    _fieldAvg->allocate();

    _fieldAvg->label(fieldIn.label());
    switch (fieldIn.vectorFieldType())
      { // switch
      case FieldBase::MULTI_SCALAR:
	_fieldAvg->vectorFieldType(FieldBase::SCALAR);
	break;
      case FieldBase::MULTI_VECTOR:
	_fieldAvg->vectorFieldType(FieldBase::VECTOR);
	break;
      case FieldBase::MULTI_TENSOR:
	_fieldAvg->vectorFieldType(FieldBase::TENSOR);
	break;
      case FieldBase::MULTI_OTHER:
	_fieldAvg->vectorFieldType(FieldBase::OTHER);
	break;
      case FieldBase::SCALAR:
      case FieldBase::VECTOR:
      case FieldBase::TENSOR:
      case FieldBase::OTHER:
      default :
	std::cerr << "Bad vector field type for CellFilterAvg." << std::endl;
	assert(0);
      } // switch
  } // if
  assert(0 != _fieldAvg);
  const ALE::Obj<RealSection>& sectionAvg = _fieldAvg->section();

  double_array fieldAvgCell(fiberDim);
  double scalar = 0.0;
  for (int iQuad=0; iQuad < numQuadPts; ++iQuad)
    scalar += wts[iQuad];

  // Loop over cells
  for (SieveMesh::label_sequence::iterator c_iter=cells->begin();
       c_iter != cellsEnd;
       ++c_iter) {
    const double* values = sectionIn->restrictPoint(*c_iter);
    
    fieldAvgCell = 0.0;
    for (int iQuad=0; iQuad < numQuadPts; ++iQuad)
      for (int i=0; i < fiberDim; ++i)
	fieldAvgCell[i] += wts[iQuad] / scalar * values[iQuad*fiberDim+i];

    _sectionAvg->updatePoint(*c_iter, &fieldAvgCell[0]);
    PetscLogFlops( numQuadPts*fiberDim*3 );
  } // for

  return *_fieldAvg;
} // filter


// End of file
