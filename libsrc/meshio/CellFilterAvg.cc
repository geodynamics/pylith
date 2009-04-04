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
pylith::meshio::CellFilterAvg<mesh_type>::CellFilterAvg(void) :
  _fieldAvg(0)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor
template<typename mesh_type>
pylith::meshio::CellFilterAvg<mesh_type>::~CellFilterAvg(void)
{ // destructor
  delete _fieldAvg; _fieldAvg = 0;
} // destructor  

// ----------------------------------------------------------------------
// Copy constructor.
template<typename mesh_type>
pylith::meshio::CellFilterAvg<mesh_type>::CellFilterAvg(const CellFilterAvg& f) :
  CellFilter<mesh_type>(f),
  _fieldAvg(0)
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
  typedef typename mesh_type::RealSection RealSection;
  typedef typename mesh_type::SieveMesh SieveMesh;
  typedef typename SieveMesh::label_sequence label_sequence;

  assert(0 != CellFilter<mesh_type>::_quadrature);

  const int numQuadPts = CellFilter<mesh_type>::_quadrature->numQuadPts();
  const double_array& wts = CellFilter<mesh_type>::_quadrature->quadWts();
  
  const ALE::Obj<SieveMesh>& sieveMesh = fieldIn.mesh().sieveMesh();
  assert(!sieveMesh.isNull());

  const ALE::Obj<label_sequence>& cells = (0 == label) ?
    sieveMesh->heightStratum(0) :
    sieveMesh->getLabelStratum(label, labelId);
  assert(!cells.isNull());
  const typename label_sequence::iterator cellsEnd = cells->end();

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
    _fieldAvg->newSection(sectionIn->getChart(), fiberDim);
    _fieldAvg->allocate();

    _fieldAvg->label(fieldIn.label());
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
  for (typename label_sequence::iterator c_iter=cells->begin();
       c_iter != cellsEnd;
       ++c_iter) {
    const double* values = sectionIn->restrictPoint(*c_iter);
    
    fieldAvgCell = 0.0;
    for (int iQuad=0; iQuad < numQuadPts; ++iQuad)
      for (int i=0; i < fiberDim; ++i)
	fieldAvgCell[i] += wts[iQuad] / scalar * values[iQuad*fiberDim+i];

    sectionAvg->updatePoint(*c_iter, &fieldAvgCell[0]);
    PetscLogFlops( numQuadPts*fiberDim*3 );
  } // for

  return *_fieldAvg;
} // filter


// End of file
