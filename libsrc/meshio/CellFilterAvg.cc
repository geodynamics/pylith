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

#include "CellFilterAvg.hh" // implementation of class methods

#include "pylith/feassemble/Quadrature.hh" // USES Quadrature

// ----------------------------------------------------------------------
// Constructor
pylith::meshio::CellFilterAvg::CellFilterAvg(void)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor
pylith::meshio::CellFilterAvg::~CellFilterAvg(void)
{ // destructor
} // destructor  

// ----------------------------------------------------------------------
// Copy constructor.
pylith::meshio::CellFilterAvg::CellFilterAvg(const CellFilterAvg& f) :
  CellFilter(f)
{ // copy constructor
} // copy constructor

// ----------------------------------------------------------------------
// Create copy of filter.
pylith::meshio::CellFilter*
pylith::meshio::CellFilterAvg::clone(void) const
{ // clone
  return new CellFilterAvg(*this);
} // clone

// ----------------------------------------------------------------------
// Filter field.
const ALE::Obj<pylith::real_section_type>&
pylith::meshio::CellFilterAvg::filter(
				  VectorFieldEnum* fieldType,
				  const ALE::Obj<real_section_type>& fieldIn,
				  const ALE::Obj<Mesh>& mesh,
				  const char* label,
				  const int labelId)
{ // filter
  assert(0 != fieldType);
  assert(0 != _quadrature);

  const int numQuadPts = _quadrature->numQuadPts();
  const double_array& wts = _quadrature->quadWts();
  
  const ALE::Obj<Mesh::label_sequence>& cells = (0 == label) ?
    mesh->heightStratum(0) :
    mesh->getLabelStratum(label, labelId);
  assert(!cells.isNull());
  const Mesh::label_sequence::iterator cellsEnd = cells->end();

  // Only processors with cells for output get the correct fiber dimension.
  const int totalFiberDim = fieldIn->getFiberDimension(*cells->begin());
  const int fiberDim = totalFiberDim / numQuadPts;
  assert(fiberDim * numQuadPts == totalFiberDim);

  *fieldType = OTHER_FIELD; // Don't know field type
  
  // Allocation field if necessary
  if (_fieldAvg.isNull() ||
      fiberDim != _fieldAvg->getFiberDimension(*cells->begin())) {
    _fieldAvg = new real_section_type(mesh->comm(), mesh->debug());
    _fieldAvg->setChart(real_section_type::chart_type(*std::min_element(cells->begin(), cells->end()),
                                                      *std::max_element(cells->begin(), cells->end())+1));
    _fieldAvg->setFiberDimension(cells, fiberDim);
    mesh->allocate(_fieldAvg);
  } // if

  double_array fieldAvgCell(fiberDim);
  double scalar = 0.0;
  for (int iQuad=0; iQuad < numQuadPts; ++iQuad)
    scalar += wts[iQuad];

  // Loop over cells
  for (Mesh::label_sequence::iterator c_iter=cells->begin();
       c_iter != cellsEnd;
       ++c_iter) {
    const real_section_type::value_type* values = 
      fieldIn->restrictPoint(*c_iter);
    
    fieldAvgCell = 0.0;
    for (int iQuad=0; iQuad < numQuadPts; ++iQuad)
      for (int i=0; i < fiberDim; ++i)
	fieldAvgCell[i] += wts[iQuad] / scalar * values[iQuad*fiberDim+i];

    _fieldAvg->updatePoint(*c_iter, &fieldAvgCell[0]);
    PetscLogFlops( numQuadPts*fiberDim*3 );
  } // for

  return _fieldAvg;
} // filter


// End of file
