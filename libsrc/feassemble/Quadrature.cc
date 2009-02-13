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

//#include "Quadrature.hh" // implementation of class methods

#include "CellGeometry.hh" // USES CellGeometry

#if 0
#include "Quadrature0D.hh"
#include "Quadrature1D.hh"
#include "Quadrature1Din2D.hh"
#include "Quadrature1Din3D.hh"
#include "Quadrature2D.hh"
#include "Quadrature2Din3D.hh"
#include "Quadrature3D.hh"
#endif

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
  _quadPtsField(0),
  _jacobianField(0),
  _jacobianDetField(0),
  _basisDerivField(0),
  _checkConditioning(false)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor
template<typename mesh_type>
pylith::feassemble::Quadrature<mesh_type>::~Quadrature(void)
{ // destructor
  delete _engine; _engine = 0;
  delete _quadPtsField; _quadPtsField = 0;
  delete _jacobianField; _jacobianField = 0;
  delete _jacobianDetField; _jacobianDetField = 0;
  delete _basisDerivField; _basisDerivField = 0;
} // destructor
  
// ----------------------------------------------------------------------
// Copy constructor
template<typename mesh_type>
pylith::feassemble::Quadrature<mesh_type>::Quadrature(const Quadrature& q) :
  QuadratureBase(q),
  _engine(0),
  _quadPtsField(0),
  _jacobianField(0),
  _jacobianDetField(0),
  _basisDerivField(0),
  _checkConditioning(q._checkConditioning)
{ // copy constructor
  if (0 != q._engine)
    _engine = q._engine->clone();
} // copy constructor

// ----------------------------------------------------------------------
template<typename mesh_type>
void
pylith::feassemble::Quadrature<mesh_type>::computeGeometry(
       const mesh_type& mesh,
       const ALE::Obj<typename mesh_type::SieveMesh::label_sequence>& cells)
{ // precomputeGeometry
  typedef typename mesh_type::RealSection RealSection;

  _setupEngine();

  const char* loggingStage = "QuadratureCreation";
  ALE::MemoryLogger& logger = ALE::MemoryLogger::singleton();
  logger.stagePush(loggingStage);

  clear();

  // Allocate field and cell buffer for quadrature points
  int fiberDim = _numQuadPts * _spaceDim;
  _quadPtsField = new topology::Field<mesh_type>(mesh);
  assert(0 != _quadPtsField);
  _quadPtsField->newSection(cells, fiberDim);
  _quadPtsField->allocate();

  // Get chart for reuse in other fields
  const ALE::Obj<RealSection>& section = _quadPtsField->section(); 
  assert(!section.isNull());
  const typename RealSection::chart_type& chart = section->getChart();

  // Allocate field and cell buffer for Jacobian at quadrature points
  fiberDim = _numQuadPts * _cellDim * _spaceDim;
  _jacobianField = new topology::Field<mesh_type>(mesh);
  assert(0 != _jacobianField);
  _jacobianField->newSection(chart, fiberDim);
  _jacobianField->allocate();
  
  // Allocate field and cell buffer for determinant of Jacobian at quad pts
  fiberDim = _numQuadPts;
  _jacobianDetField = new topology::Field<mesh_type>(mesh);
  assert(0 != _jacobianDetField);
  _jacobianDetField->newSection(chart, fiberDim);
  _jacobianDetField->allocate();
  
  // Allocate field for derivatives of basis functions at quad pts
  fiberDim = _numQuadPts * _numBasis * _spaceDim;
  _basisDerivField = new topology::Field<mesh_type>(mesh);
  assert(0 != _basisDerivField);
  _basisDerivField->newSection(chart, fiberDim);
  _basisDerivField->allocate();

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

  typedef typename mesh_type::SieveMesh::label_sequence label_sequence;
  typedef ALE::ISieveVisitor::RestrictVisitor<RealSection> RealSectionVisitor;

  const typename label_sequence::iterator cellsEnd = cells->end();
  assert(0 != _geometry);
  const int numCorners = _geometry->numCorners();
  const ALE::Obj<typename mesh_type::SieveMesh>& sieveMesh = mesh.sieveMesh();
  assert(!sieveMesh.isNull());
  const ALE::Obj<RealSection>& coordinates = 
    sieveMesh->getRealSection("coordinates");
  RealSectionVisitor coordsVisitor(coordinates, numCorners*_spaceDim);

  const ALE::Obj<RealSection>& quadPtsSection = _quadPtsField->section();
  const ALE::Obj<RealSection>& jacobianSection = _jacobianField->section();
  const ALE::Obj<RealSection>& jacobianDetSection = 
    _jacobianDetField->section();
  const ALE::Obj<RealSection>& basisDerivSection = _basisDerivField->section();

  for(typename label_sequence::iterator c_iter = cells->begin();
      c_iter != cellsEnd;
      ++c_iter) {
    sieveMesh->restrictClosure(*c_iter, coordsVisitor);
    const double* cellVertexCoords = coordsVisitor.getValues();
    assert(0 != cellVertexCoords);
    _resetGeometry();
    computeGeometry(cellVertexCoords, _spaceDim, *c_iter);

    // Update fields with cell data
    quadPtsSection->updatePoint(*c_iter, &_quadPts[0]);
    jacobianSection->updatePoint(*c_iter, &_jacobian[0]);
    jacobianDetSection->updatePoint(*c_iter, &_jacobianDet[0]);
    basisDerivSection->updatePoint(*c_iter, &_basisDeriv[0]);
  } // for
} // computeGeometry

// ----------------------------------------------------------------------
template<typename mesh_type>
void
pylith::feassemble::Quadrature<mesh_type>::retrieveGeometry(const typename mesh_type::SieveMesh::point_type& cell)
{ // retrieveGeometry
  typedef typename mesh_type::RealSection RealSection;

  assert(0 != _quadPtsField);
  assert(0 != _jacobianField);
  assert(0 != _jacobianDetField);
  assert(0 != _basisDerivField);

  const ALE::Obj<RealSection>& quadPtsSection = _quadPtsField->section();
  quadPtsSection->restrictPoint(cell, &_quadPts[0], _quadPts.size());

  const ALE::Obj<RealSection>& jacobianSection = _jacobianField->section();
  jacobianSection->restrictPoint(cell, &_jacobian[0], _jacobian.size());

  const ALE::Obj<RealSection>& jacobianDetSection = 
    _jacobianDetField->section();
  jacobianDetSection->restrictPoint(cell, 
				    &_jacobianDet[0], _jacobianDet.size());

  const ALE::Obj<RealSection>& basisDerivSection = _basisDerivField->section();
  basisDerivSection->restrictPoint(cell, &_basisDeriv[0], _basisDeriv.size());
} // retrieveGeometry

// ----------------------------------------------------------------------
// Deallocate temporary storage;
template<typename mesh_type>
void
pylith::feassemble::Quadrature<mesh_type>::clear(void)
{ // clear
  delete _engine; _engine = 0;

  // Clear storage for fields
  delete _quadPtsField; _quadPtsField = 0;
  delete _jacobianField; _jacobianField = 0;
  delete _jacobianDetField; _jacobianDetField = 0;
  delete _basisDerivField; _basisDerivField = 0;
} // clear

// ----------------------------------------------------------------------
// Setup quadrature engine.
template<typename mesh_type>
void
pylith::feassemble::Quadrature<mesh_type>::_setupEngine(void)
{ // clear
  delete _engine; _engine = 0;

  const int cellDim = _cellDim;
  const int spaceDim = _spaceDim;

#if 0
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
  else if (2 == spaceDim) {
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
  else if (3 == spaceDim) {
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
#endif
} // _setupEngine


// End of file 
