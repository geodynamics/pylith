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

#include "CellGeometry.hh" // USES CellGeometry

#include "QuadratureEngine.hh" // USES QuadratureEngine
#include "Quadrature0D.hh"
#include "Quadrature1D.hh"
#include "Quadrature1Din2D.hh"
#include "Quadrature1Din3D.hh"
#include "Quadrature2D.hh"
#include "Quadrature2Din3D.hh"
#include "Quadrature3D.hh"

#include "pylith/topology/Fields.hh" // HOLDSA Fields
#include "pylith/topology/Field.hh" // HOLDSA Field

#include <cstring> // USES memcpy()
#include <cassert> // USES assert()
#include <stdexcept> // USES std::runtime_error
#include <sstream> // USES std::ostringstream

// ----------------------------------------------------------------------
// Constructor
template<typename mesh_type>
pylith::feassemble::Quadrature<mesh_type>::Quadrature(void) :
  _engine(0),
  _geometryFields(0),
  _checkConditioning(false)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor
template<typename mesh_type>
pylith::feassemble::Quadrature<mesh_type>::~Quadrature(void)
{ // destructor
  deallocate();
} // destructor
  
// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
template<typename mesh_type>
void
pylith::feassemble::Quadrature<mesh_type>::deallocate(void)
{ // deallocate
  QuadratureRefCell::deallocate();

  delete _engine; _engine = 0;
  delete _geometryFields; _geometryFields = 0;
} // deallocate
  
// ----------------------------------------------------------------------
// Copy constructor
template<typename mesh_type>
pylith::feassemble::Quadrature<mesh_type>::Quadrature(const Quadrature& q) :
  QuadratureRefCell(q),
  _engine(0),
  _geometryFields(0),
  _checkConditioning(q._checkConditioning)
{ // copy constructor
  if (0 != q._engine)
    _engine = q._engine->clone();
} // copy constructor

// ----------------------------------------------------------------------
// Setup quadrature engine.
template<typename mesh_type>
void
pylith::feassemble::Quadrature<mesh_type>::initializeGeometry(void)
{ // initializeGeometry
  clear();
  assert(0 == _engine);

  const int cellDim = _cellDim;
  const int spaceDim = _spaceDim;

  if (1 == spaceDim)
    if (1 == cellDim)
      _engine = new Quadrature1D(*this);
    else if (0 == cellDim)
      _engine = new Quadrature0D(*this);
    else {
      std::cerr << "Unknown quadrature case with cellDim '" 
		<< cellDim << "' and spaceDim '" << spaceDim << "'" 
		<< std::endl;
      assert(0);
    } // if/else
  else if (2 == spaceDim)
    if (2 == cellDim)
      _engine = new Quadrature2D(*this);
    else if (1 == cellDim)
      _engine = new Quadrature1Din2D(*this);
    else if (0 == cellDim)
      _engine = new Quadrature0D(*this);
    else {
      std::cerr << "Unknown quadrature case with cellDim '" 
		<< cellDim << "' and spaceDim '" << spaceDim << "'" 
		<< std::endl;
      assert(0);
    } // if/else
  else if (3 == spaceDim)
    if (3 == cellDim)
      _engine = new Quadrature3D(*this);
    else if (2 == cellDim)
      _engine = new Quadrature2Din3D(*this);
    else if (1 == cellDim)
      _engine = new Quadrature1Din3D(*this);
    else if (0 == cellDim)
      _engine = new Quadrature0D(*this);
    else {
      std::cerr << "Unknown quadrature case with cellDim '" 
		<< cellDim << "' and spaceDim '" << spaceDim << "'" 
		<< std::endl;
      assert(0);
    } // if/else
  else {
    std::cerr << "Unknown quadrature case with cellDim '" 
	      << cellDim << "' and spaceDim '" << spaceDim << "'" 
	      << std::endl;
    assert(0);
  } // if/else

  assert(0 != _engine);
  _engine->initialize();
} // initializeGeometry

// ----------------------------------------------------------------------
// Compute geometric quantities for each cell.
template<typename mesh_type>
void
pylith::feassemble::Quadrature<mesh_type>::computeGeometry(
       const mesh_type& mesh,
       const ALE::Obj<typename mesh_type::SieveMesh::label_sequence>& cells)
{ // computeGeometry
  assert(0 != _engine);

  typedef typename mesh_type::RealSection RealSection;
  typedef typename mesh_type::SieveMesh::label_sequence label_sequence;
  typedef typename mesh_type::RestrictVisitor RestrictVisitor;

  const char* loggingStage = "Quadrature";
  ALE::MemoryLogger& logger = ALE::MemoryLogger::singleton();
  logger.stagePush(loggingStage);

  delete _geometryFields;
  _geometryFields = new topology::Fields<topology::Field<mesh_type> >;

  // Allocate field and cell buffer for quadrature points
  _geometryFields->add("quadrature points", "quadrature_points");
  topology::Field<mesh_type>& quadPtsField = 
    _geometryFields->get("quadrature points");
  int fiberDim = _numQuadPts * _spaceDim;
  quadPtsField.newSection(cells, fiberDim);
  quadPtsField.allocate();

  // Get chart for reuse in other fields
  const ALE::Obj<RealSection>& section = quadPtsField->section(); 
  assert(!section.isNull());
  const typename RealSection::chart_type& chart = section->getChart();

  // Allocate field and cell buffer for Jacobian at quadrature points
  std::cout << "Jacobian: cell dim: " << _cellDim << std::endl;
  _geometryFields->add("jacobian", "jacobian");
  topology::Field<mesh_type>& jacobianField = 
    _geometryFields->get("jacobian");
  fiberDim = (_cellDim > 0) ?
    _numQuadPts * _cellDim * _spaceDim :
    _numQuadPts * 1 * _spaceDim;
  jacobianField->newSection(chart, fiberDim);
  jacobianField->allocate();
  
  // Allocate field and cell buffer for determinant of Jacobian at quad pts
  std::cout << "Jacobian det:" << std::endl;
  _geometryFields->add("determinant(jacobian)", "determinant_jacobian");
  topology::Field<mesh_type>& jacobianDetField = 
    _geometryFields->get("determinant(jacobian)");
  fiberDim = _numQuadPts;
  jacobianDetField.newSection(chart, fiberDim);
  jacobianDetField.allocate();
  
  // Allocate field for derivatives of basis functions at quad pts
  std::cout << "Basis derivatives: num basis: " << _numBasis << std::endl;
  _geometryFields->add("derivative basis functions",
		       "derivative_basis_functions");
  topology::Field<mesh_type>& basisDerivField = 
    _geometryFields->get("jacobian");
  fiberDim = _numQuadPts * _numBasis * _spaceDim;
  basisDerivField.newSection(chart, fiberDim);
  basisDerivField.allocate();

  logger.stagePop();

#if defined(ALE_MEM_LOGGING)
  std::cout 
    << loggingStage << ": " 
    << logger.getNumAllocations(loggingStage)
    << " allocations " << logger.getAllocationTotal(loggingStage)
    << " bytes"
    << std::endl
    << loggingStage << ": "
    << logger.getNumDeallocations(loggingStage)
    << " deallocations " << logger.getDeallocationTotal(loggingStage)
    << " bytes"
    << std::endl;
#endif

  assert(!cells.isNull());
  const typename label_sequence::iterator cellsBegin = cells->begin();
  const typename label_sequence::iterator cellsEnd = cells->end();
  assert(0 != _geometry);
  const int numBasis = _numBasis;
  const ALE::Obj<typename mesh_type::SieveMesh>& sieveMesh = mesh.sieveMesh();
  assert(!sieveMesh.isNull());

  scalar_array coordinatesCell(numBasis*_spaceDim);
  const ALE::Obj<RealSection>& coordinates = 
    sieveMesh->getRealSection("coordinates");
  RestrictVisitor coordsVisitor(*coordinates,
				coordinatesCell.size(), &coordinatesCell[0]);

  const ALE::Obj<RealSection>& quadPtsSection = quadPtsField->section();
  const ALE::Obj<RealSection>& jacobianSection = jacobianField->section();
  const ALE::Obj<RealSection>& jacobianDetSection = 
    jacobianDetField->section();
  const ALE::Obj<RealSection>& basisDerivSection = basisDerivField->section();

  const scalar_array& quadPts = _engine->quadPts();
  const scalar_array& jacobian = _engine->jacobian();
  const scalar_array& jacobianDet = _engine->jacobianDet();
  const scalar_array& basisDeriv = _engine->basisDeriv();

  for(typename label_sequence::iterator c_iter = cellsBegin;
      c_iter != cellsEnd;
      ++c_iter) {
    coordsVisitor.clear();
    sieveMesh->restrictClosure(*c_iter, coordsVisitor);
    _engine->computeGeometry(coordinatesCell, *c_iter);

    // Update fields with cell data
    quadPtsSection->updatePoint(*c_iter, &quadPts[0]);
    jacobianSection->updatePoint(*c_iter, &jacobian[0]);
    jacobianDetSection->updatePoint(*c_iter, &jacobianDet[0]);
    basisDerivSection->updatePoint(*c_iter, &basisDeriv[0]);
  } // for
} // computeGeometry

// ----------------------------------------------------------------------
template<typename mesh_type>
void
pylith::feassemble::Quadrature<mesh_type>::retrieveGeometry(const typename mesh_type::SieveMesh::point_type& cell)
{ // retrieveGeometry
  assert(0 != _geometryFields);
  assert(0 != _engine);

  typedef typename mesh_type::RealSection RealSection;

  const scalar_array& quadPts = _engine->quadPts();
  const scalar_array& jacobian = _engine->jacobian();
  const scalar_array& jacobianDet = _engine->jacobianDet();
  const scalar_array& basisDeriv = _engine->basisDeriv();

  const ALE::Obj<RealSection>& quadPtsSection =
    _geometryFields->get("quadrature points").section();
  quadPtsSection->restrictPoint(cell, const_cast<PylithScalar*>(&quadPts[0]),
				quadPts.size());

  const ALE::Obj<RealSection>& jacobianSection =
    _geometryFields->get("jacobian").section();
  jacobianSection->restrictPoint(cell, const_cast<PylithScalar*>(&jacobian[0]),
				 jacobian.size());

  const ALE::Obj<RealSection>& jacobianDetSection = 
    _geometryFields->get("determinant(jacobian)").section();
  jacobianDetSection->restrictPoint(cell, const_cast<PylithScalar*>(&jacobianDet[0]),
				    jacobianDet.size());

  const ALE::Obj<RealSection>& basisDerivSection =
    _geometryFields->get("determinant basisfunctions").section();
  basisDerivSection->restrictPoint(cell, const_cast<PylithScalar*>(&basisDeriv[0]),
				   basisDeriv.size());
} // retrieveGeometry

// ----------------------------------------------------------------------
// Deallocate temporary storage;
template<typename mesh_type>
void
pylith::feassemble::Quadrature<mesh_type>::clear(void)
{ // clear
  delete _engine; _engine = 0;

  // Clear storage of precomputed geometry.
  delete _geometryFields; _geometryFields = 0;
} // clear


// End of file 
