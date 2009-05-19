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

#include "pylith/topology/Field.hh" // USES Field

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
  delete _fieldVecNorm; _fieldVecNorm = 0;
} // destructor  

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
const field_type&
pylith::meshio::VertexFilterVecNorm<field_type>::filter(
				   const field_type& fieldIn)
{ // filter
  typedef typename field_type::Mesh::RealSection RealSection;
  typedef typename field_type::Mesh::SieveMesh SieveMesh;
  typedef typename SieveMesh::label_sequence label_sequence;

  const ALE::Obj<SieveMesh>& sieveMesh = fieldIn.mesh().sieveMesh();
  assert(!sieveMesh.isNull());

  const ALE::Obj<label_sequence>& vertices = sieveMesh->depthStratum(0);
  assert(!vertices.isNull());
  const typename label_sequence::iterator verticesEnd = vertices->end();

  const ALE::Obj<RealSection>& sectionIn = fieldIn.section();
  assert(!sectionIn.isNull());
  const int fiberDimIn = (vertices->size() > 0) ? 
    sectionIn->getFiberDimension(*vertices->begin()) : 0;
  const int fiberDimNorm = 1;

  // Allocation field if necessary
  if (0 == _fieldVecNorm) {
    _fieldVecNorm = new field_type(fieldIn.mesh());
    _fieldVecNorm->label("vector norm");
    _fieldVecNorm->newSection(sectionIn->getChart(), fiberDimNorm);
    _fieldVecNorm->allocate();

    _fieldVecNorm->label(fieldIn.label());    
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
	std::cerr << "Bad vector field type for VertexFilterVecNorm." << std::endl;
	assert(0);
      } // switch
  } // if

  const ALE::Obj<RealSection>& sectionNorm = 
    _fieldVecNorm->section();
  assert(!sectionNorm.isNull());

  double norm = 0.0;
  // Loop over vertices
  for (typename label_sequence::iterator v_iter=vertices->begin();
       v_iter != verticesEnd;
       ++v_iter) {
    const double* values = sectionIn->restrictPoint(*v_iter);
    
    norm = 0.0;
    for (int i=0; i < fiberDimIn; ++i)
      norm += values[i]*values[i];
    norm = sqrt(norm);

    sectionNorm->updatePoint(*v_iter, &norm);
  } // for
  PetscLogFlops(vertices->size() * (1 + 2*fiberDimIn) );

  return *_fieldVecNorm;
} // filter


// End of file
