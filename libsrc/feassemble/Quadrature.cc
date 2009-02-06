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

#include "Quadrature.hh" // implementation of class methods

#include "CellGeometry.hh" // USES CellGeometry

#include <cstring> // USES memcpy()
#include <cassert> // USES assert()
#include <stdexcept> // USES std::runtime_error
#include <sstream> // USES std::ostringstream

// ----------------------------------------------------------------------
// Constructor
pylith::feassemble::Quadrature::Quadrature(void) :
  _minJacobian(0),
  _cellDim(0),
  _numBasis(0),
  _numQuadPts(0),
  _spaceDim(0),
  _geometry(0),
  _precomputed(false),
  _checkConditioning(false)
{ // constructor
  _quadPtsPre     = new real_section_type(PETSC_COMM_WORLD);
  _jacobianPre    = new real_section_type(PETSC_COMM_WORLD);
  _jacobianDetPre = new real_section_type(PETSC_COMM_WORLD);
  _jacobianInvPre = new real_section_type(PETSC_COMM_WORLD);
  _basisDerivPre  = new real_section_type(PETSC_COMM_WORLD);
} // constructor

// ----------------------------------------------------------------------
// Destructor
pylith::feassemble::Quadrature::~Quadrature(void)
{ // destructor
  delete _geometry; _geometry = 0;
} // destructor
  
// ----------------------------------------------------------------------
// Copy constructor
pylith::feassemble::Quadrature::Quadrature(const Quadrature& q) :
  _minJacobian(q._minJacobian),
  _quadPtsRef(q._quadPtsRef),
  _quadPts(q._quadPts),
  _quadWts(q._quadWts),
  _basis(q._basis),
  _basisDerivRef(q._basisDerivRef),
  _basisDeriv(q._basisDeriv),
  _jacobian(q._jacobian),
  _jacobianInv(q._jacobianInv),
  _jacobianDet(q._jacobianDet),
  _cellDim(q._cellDim),
  _numBasis(q._numBasis),
  _numQuadPts(q._numQuadPts),
  _spaceDim(q._spaceDim),
  _geometry(0),
  _precomputed(q._precomputed),
  _checkConditioning(q._checkConditioning)
{ // copy constructor
  if (0 != q._geometry)
    _geometry = q._geometry->clone();
  _quadPtsPre = q._quadPtsPre;
  _jacobianPre = q._jacobianPre;
  _jacobianDetPre = q._jacobianDetPre;
  _jacobianInvPre = q._jacobianInvPre;
  _basisDerivPre = q._basisDerivPre;
} // copy constructor

// ----------------------------------------------------------------------
// Set basis functions and their derivatives and coordinates and
//   weights of the quadrature points.
void
pylith::feassemble::Quadrature::initialize(const double* basis,
					   const double* basisDerivRef,
					   const double* quadPtsRef,
					   const double* quadWts,
					   const int cellDim,
					   const int numBasis,
					   const int numQuadPts,
					   const int spaceDim)
{ // initialize
  if (0 == basis ||
      0 == basisDerivRef ||
      0 == quadPtsRef ||
      0 == quadWts ||
      cellDim < 0 || cellDim > 3 ||
      numBasis < 1 ||
      numQuadPts < 1 ||
      spaceDim < 1 || spaceDim > 3) {
    std::ostringstream msg;
    msg << "Incompatible values for quadrature information. Basis functions,\n"
	<< "their derivatives, and coordinates and weights of quadrature\n"
	<< "points must all be specified.\n"
	<< "Values:\n"
	<< "  basis pointer: " << basis << "\n"
	<< "  basis derivatites pointer: " << basisDerivRef << "\n"
	<< "  quadrature points pointer: " << quadPtsRef << "\n"
	<< "  quadrature weights pointer: " << quadWts << "\n"
	<< "  space dimension: " << spaceDim << "\n"
	<< "  # basis functions: " << numBasis << "\n"
	<< "  # quadrature points: " << numQuadPts << "\n"
	<< "  dimension of coordinate space: " << spaceDim << "\n";
    throw std::runtime_error(msg.str());
  } // if

  if (cellDim > 0) {
    int size = numBasis * numQuadPts; assert(size > 0);
    _basis.resize(size);
    for (int i=0; i < size; ++i)
      _basis[i] = basis[i];

    size = numQuadPts * numBasis * cellDim; assert(size > 0);
    _basisDerivRef.resize(size);
    for (int i=0; i < size; ++i)
      _basisDerivRef[i] = basisDerivRef[i];

    size = numQuadPts * cellDim; assert(size > 0);
    _quadPtsRef.resize(size);
    for (int i=0; i < size; ++i)
      _quadPtsRef[i] = quadPtsRef[i];

    size = numQuadPts; assert(size > 0);
    _quadWts.resize(size);
    for (int i=0; i < size; ++i)
      _quadWts[i] = quadWts[i];

    _cellDim = cellDim;
    _numBasis = numBasis;
    _numQuadPts = numQuadPts;
    _spaceDim = spaceDim;

    // Allocate for Jacobian and its inverse
    size = numQuadPts * cellDim * spaceDim; assert(size > 0);
    _jacobian.resize(size);
    _jacobianInv.resize(size);

    // Allocate for Jacobian determinant
    size = numQuadPts; assert(size > 0);
    _jacobianDet.resize(size);

    // Allocate for basis derivatives (in global coordinates)
    size = numQuadPts * numBasis * spaceDim; assert(size > 0);
    _basisDeriv.resize(size);

    // Allocate for quad pts
    size = numQuadPts*spaceDim; assert(size > 0);
    _quadPts.resize(size);
  } else {
    if (1 != numBasis ||
	1 != numQuadPts ||
	1 != spaceDim) {
      std::ostringstream msg;
      msg << "0-D quadrature only works in 1-D and is limited to 1 basis "
	  << "function and 1 quadrature point.\n"
	  << "Values:\n"
	  << "  cell dimension: " << cellDim << "\n"
	  << "  spatial dimension: " << spaceDim << "\n"
	  << "  # basis functions: " << numBasis << "\n"
	  << "  # quadrature points: " << numQuadPts << "\n";
      throw std::runtime_error(msg.str());
    } // if

    int size = 1;
    _basis.resize(size);
    for (int i=0; i < size; ++i)
      _basis[i] = basis[i];

    size = 1;
    _basisDerivRef.resize(size);
    for (int i=0; i < size; ++i)
      _basisDerivRef[i] = basisDerivRef[i];

    size = 1;
    _quadPtsRef.resize(size);
    for (int i=0; i < size; ++i)
      _quadPtsRef[i] = quadPtsRef[i];

    size = 1;
    _quadWts.resize(size);
    for (int i=0; i < size; ++i)
      _quadWts[i] = quadWts[i];

    _cellDim = cellDim;
    _numBasis = numBasis;
    _numQuadPts = numQuadPts;
    _spaceDim = spaceDim;

    // Allocate for Jacobian and its inverse
    size = 1;
    _jacobian.resize(size);
    _jacobianInv.resize(size);

    // Allocate for Jacobian determinant
    size = 1;
    _jacobianDet.resize(size);

    // Allocate for basis derivatives (in global coordinates)
    size = numQuadPts * numBasis * spaceDim; assert(size > 0);
    _basisDeriv.resize(size);

    // Allocate for quad pts
    size = spaceDim; assert(size > 0);
    _quadPts.resize(size);
  } // else
} // initialize

// ----------------------------------------------------------------------
// Set geometry associated with reference cell.
void
pylith::feassemble::Quadrature::refGeometry(CellGeometry* const geometry)
{ // refGeometry
  delete _geometry; _geometry = (0 != geometry) ? geometry->clone() : 0;
} // refGeometry

// ----------------------------------------------------------------------
// Get geometry associated with reference cell.
const pylith::feassemble::CellGeometry&
pylith::feassemble::Quadrature::refGeometry(void) const
{ // refGeometry
  assert(0 != _geometry);
  return *_geometry;
} // refGeometry

// ----------------------------------------------------------------------
// Set entries in geometry arrays to zero.
void
pylith::feassemble::Quadrature::_resetGeometry(void)
{ // _resetGeometry
  _jacobian = 0.0;
  _jacobianDet = 0.0;
  _jacobianInv = 0.0;
  _basisDeriv = 0.0;
  _quadPts = 0.0;
} // _resetGeometry

// ----------------------------------------------------------------------
// Check determinant of Jacobian against minimum allowable value
void
pylith::feassemble::Quadrature::_checkJacobianDet(const double det,
						  const int cell) const
{ // _checkJacobianDet
  if (det < _minJacobian) {
    std::ostringstream msg;
    msg << "Determinant of Jacobian (" << det << ") for cell " << cell
	<< " is smaller than minimum permissible value (" << _minJacobian
	<< ")!\n";
    throw std::runtime_error(msg.str());
  } // if
} // _checkJacobianDet

// ----------------------------------------------------------------------
void
pylith::feassemble::Quadrature::resetPrecomputation()
{ // resetPrecomputation
  _precomputed = false;
  _quadPtsPre->clear();
  _jacobianPre->clear();
  _jacobianDetPre->clear();
  _jacobianInvPre->clear();
  _basisDerivPre->clear();
} // resetPrecomputation

// ----------------------------------------------------------------------
void
pylith::feassemble::Quadrature::precomputeGeometry(
			const ALE::Obj<Mesh>& mesh,
			const ALE::Obj<real_section_type>& coordinates,
			const ALE::Obj<Mesh::label_sequence>& cells)
{ // precomputeGeometry
  if (_precomputed)
    return;

  _quadPtsPre->setChart(real_section_type::chart_type(
		      *std::min_element(cells->begin(), cells->end()),
		      *std::max_element(cells->begin(), cells->end())+1));
  _quadPtsPre->setFiberDimension(cells, _numQuadPts*_spaceDim);
  _quadPtsPre->allocatePoint();
  _quadPtsPreV = new ALE::ISieveVisitor::RestrictVisitor<real_section_type>(
			  *_quadPtsPre, _numQuadPts*_spaceDim, &_quadPts[0]);

  _jacobianPre->getAtlas()->setAtlas(_quadPtsPre->getAtlas()->getAtlas());
  _jacobianPre->getAtlas()->allocatePoint();
  _jacobianPre->setFiberDimension(cells, _numQuadPts*_cellDim*_spaceDim);
  _jacobianPre->allocatePoint();
  _jacobianPreV = new ALE::ISieveVisitor::RestrictVisitor<real_section_type>(
		*_jacobianPre, _numQuadPts*_cellDim*_spaceDim, &_jacobian[0]);

  _jacobianDetPre->getAtlas()->setAtlas(_quadPtsPre->getAtlas()->getAtlas());
  _jacobianDetPre->getAtlas()->allocatePoint();
  _jacobianDetPre->setFiberDimension(cells, _numQuadPts);
  _jacobianDetPre->allocatePoint();
  _jacobianDetPreV = new ALE::ISieveVisitor::RestrictVisitor<real_section_type>(*_jacobianDetPre, _numQuadPts, &_jacobianDet[0]);

  _jacobianInvPre->setAtlas(_jacobianPre->getAtlas());
  _jacobianInvPre->setFiberDimension(cells, _numQuadPts*_cellDim*_spaceDim);
  _jacobianInvPre->allocatePoint();
  _jacobianInvPreV = new ALE::ISieveVisitor::RestrictVisitor<real_section_type>(*_jacobianInvPre, _numQuadPts*_cellDim*_spaceDim, &_jacobianInv[0]);
  //_jITag = _jacobianInvPre->copyCustomAtlas(_jacobianPre, _jTag);

  _basisDerivPre->getAtlas()->setAtlas(_quadPtsPre->getAtlas()->getAtlas());
  _basisDerivPre->getAtlas()->allocatePoint();
  _basisDerivPre->setFiberDimension(cells, _numQuadPts*_numBasis*_spaceDim);
  _basisDerivPre->allocatePoint();
  _basisDerivPreV = new ALE::ISieveVisitor::RestrictVisitor<real_section_type>(*_basisDerivPre, _numQuadPts*_numBasis*_spaceDim, &_basisDeriv[0]);

#if 0
  const int ncells = cells->size();
  const int nbytes = (_numQuadPts*_spaceDim + // quadPts
		      _numQuadPts*_cellDim*_spaceDim + // jacobian
		      _numQuadPts*_cellDim*_spaceDim + // jacobianInv
		      _numQuadPts + // jacobianDet
		      _numQuadPts*_numBasis*_spaceDim // basisDeriv
		      ) * ncells * sizeof(double);
  
  std::cout << "Quadrature::precomputeGeometry() allocating "
	    << nbytes/(1024*1024) << " MB."
	    << std::endl;
#endif

  const Mesh::label_sequence::iterator cellsEnd = cells->end();
  for(Mesh::label_sequence::iterator c_iter = cells->begin();
      c_iter != cellsEnd;
      ++c_iter) {
    this->computeGeometry(mesh, coordinates, *c_iter);

    // Set coordinates of quadrature points in cell
    _quadPtsPre->updatePoint(*c_iter, &_quadPts[0]);

    // Set Jacobian at quadrature points in cell
    _jacobianPre->updatePoint(*c_iter, &_jacobian[0]);

    // Set determinant of Jacobian at quadrature points in cell
    _jacobianDetPre->updatePoint(*c_iter, &_jacobianDet[0]);

    // Set inverse of Jacobian at quadrature points in cell
    _jacobianInvPre->updatePoint(*c_iter, &_jacobianInv[0]);

    // Set derivatives of basis functions with respect to global
    _basisDerivPre->updatePoint(*c_iter, &_basisDeriv[0]);
  } // for
  _precomputed = true;
} // precomputeGeometry

// ----------------------------------------------------------------------
void
pylith::feassemble::Quadrature::retrieveGeometry(const ALE::Obj<SieveMesh>& mesh,
						 const SieveMesh::point_type& cell)
{ // retrieveGeometry
  _quadPtsPreV->clear();
  mesh->restrictClosure(cell, *_quadPtsPreV);

  _jacobianPreV->clear();
  mesh->restrictClosure(cell, *_jacobianPreV);

  _jacobianDetPreV->clear();
  mesh->restrictClosure(cell, *_jacobianDetPreV);

  _jacobianInvPreV->clear();
  mesh->restrictClosure(cell, *_jacobianInvPreV);

  _basisDerivPreV->clear();
  mesh->restrictClosure(cell, *_basisDerivPreV);
} // retrieveGeometry

// End of file 
