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

#include "VertexFilterVecNorm.hh" // implementation of class methods

// ----------------------------------------------------------------------
// Constructor
pylith::meshio::VertexFilterVecNorm::VertexFilterVecNorm(void)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor
pylith::meshio::VertexFilterVecNorm::~VertexFilterVecNorm(void)
{ // destructor
} // destructor  

// ----------------------------------------------------------------------
// Copy constructor.
pylith::meshio::VertexFilterVecNorm::VertexFilterVecNorm(const VertexFilterVecNorm& f) :
  VertexFilter(f)
{ // copy constructor
} // copy constructor

// ----------------------------------------------------------------------
// Create copy of filter.
pylith::meshio::VertexFilter*
pylith::meshio::VertexFilterVecNorm::clone(void) const
{ // clone
  return new VertexFilterVecNorm(*this);
} // clone

// ----------------------------------------------------------------------
// Filter field.
const ALE::Obj<pylith::real_section_type>&
pylith::meshio::VertexFilterVecNorm::filter(
				  VectorFieldEnum* fieldType,
				  const ALE::Obj<real_section_type>& fieldIn,
				  const ALE::Obj<ALE::Mesh>& mesh)
{ // filter
  assert(0 != fieldType);
  *fieldType = SCALAR_FIELD;

  const ALE::Obj<Mesh::label_sequence>& vertices = mesh->depthStratum(0);
  assert(!vertices.isNull());
  const Mesh::label_sequence::iterator verticesEnd = vertices->end();

  const int totalFiberDim = fieldIn->getFiberDimension(*vertices->begin());
  const int fiberDim = fieldIn->getFiberDimension(*vertices->begin());

  // Allocation field if necessary
  if (_fieldVecNorm.isNull() ||
      1 != _fieldVecNorm->getFiberDimension(*vertices->begin())) {
    _fieldVecNorm = new real_section_type(mesh->comm(), mesh->debug());
    _fieldVecNorm->setFiberDimension(vertices, 1);
    mesh->allocate(_fieldVecNorm);
  } // if

  double norm = 0.0;

  // Loop over vertices
  for (Mesh::label_sequence::iterator v_iter=vertices->begin();
       v_iter != verticesEnd;
       ++v_iter) {
    const real_section_type::value_type* values = 
      fieldIn->restrictPoint(*v_iter);
    
    norm = 0.0;
    for (int i=0; i < fiberDim; ++i)
      norm += values[i]*values[i];
    norm = sqrt(norm);

    _fieldVecNorm->updatePoint(*v_iter, &norm);
  } // for

  return _fieldVecNorm;
} // filter


// End of file
