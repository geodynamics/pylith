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

#include <string.h> // USES memcpy()
#include <assert.h> // USES assert()
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
  _precomputed(false)
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
  _precomputed(q._precomputed)
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
					   const Mesh::point_type& cell) const
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
  if (_precomputed) return;
  const Mesh::label_sequence::iterator end = cells->end();

  _quadPtsPre->setFiberDimension(cells, _numQuadPts*_spaceDim);
  _quadPtsPre->allocatePoint();
  _qTag = mesh->calculateCustomAtlas(_quadPtsPre, cells);
  _jacobianPre->getAtlas()->setAtlas(_quadPtsPre->getAtlas()->getAtlas());
  _jacobianPre->setFiberDimension(cells, _numQuadPts*_cellDim*_spaceDim);
  _jacobianPre->allocatePoint();
  _jTag = mesh->calculateCustomAtlas(_jacobianPre, cells);
  _jacobianDetPre->getAtlas()->setAtlas(_quadPtsPre->getAtlas()->getAtlas());
  _jacobianDetPre->setFiberDimension(cells, _numQuadPts);
  _jacobianDetPre->allocatePoint();
  _jDTag = mesh->calculateCustomAtlas(_jacobianDetPre, cells);
  _jacobianInvPre->setAtlas(_jacobianPre->getAtlas());
  _jacobianInvPre->setFiberDimension(cells, _numQuadPts*_cellDim*_spaceDim);
  _jacobianInvPre->allocatePoint();
  _jITag = _jacobianInvPre->copyCustomAtlas(_jacobianPre, _jTag);
  _basisDerivPre->getAtlas()->setAtlas(_quadPtsPre->getAtlas()->getAtlas());
  _basisDerivPre->setFiberDimension(cells, _numQuadPts*_numBasis*_spaceDim);
  _basisDerivPre->allocatePoint();
  _bTag = mesh->calculateCustomAtlas(_basisDerivPre, cells);

  for(Mesh::label_sequence::iterator c_iter = cells->begin();
      c_iter != end;
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
pylith::feassemble::Quadrature::retrieveGeometry(
			      const ALE::Obj<Mesh>& mesh,
                  const ALE::Obj<real_section_type>& coordinates,
			      const Mesh::point_type& cell,
                  const int c)
{ // retrieveGeometry
#define FASTER
#ifdef FASTER
  const real_section_type::value_type* values =
    mesh->restrict(_quadPtsPre, _qTag, c);
#else
  const real_section_type::value_type* values =
    _quadPtsPre->restrictPoint(cell);
#endif
  int size = _numQuadPts * _spaceDim;
  assert(size == _quadPtsPre->getFiberDimension(cell));
  for(int i=0; i < size; ++i)
    _quadPts[i] = values[i];

#ifdef FASTER
  values = mesh->restrict(_jacobianPre, _jTag, c);
#else
  values = _jacobianPre->restrictPoint(cell);
#endif
  size = _numQuadPts * _cellDim * _spaceDim;
  assert(size == _jacobianPre->getFiberDimension(cell));
  for(int i=0; i < size; ++i)
    _jacobian[i] = values[i];

#ifdef FASTER
  values = mesh->restrict(_jacobianDetPre, _jDTag, c);
#else
  values = _jacobianDetPre->restrictPoint(cell);
#endif
  size = _numQuadPts;
  assert(size == _jacobianDetPre->getFiberDimension(cell));
  for(int i=0; i < size; ++i)
    _jacobianDet[i] = values[i];

#ifdef FASTER
  values = mesh->restrict(_jacobianInvPre, _jITag, c);
#else
  values = _jacobianInvPre->restrictPoint(cell);
#endif
  size = _numQuadPts * _cellDim * _spaceDim;
  assert(size == _jacobianInvPre->getFiberDimension(cell));
  for(int i=0; i < size; ++i)
    _jacobianInv[i] = values[i];

#ifdef FASTER
  values = mesh->restrict(_basisDerivPre, _bTag, c);
#else
  values = _basisDerivPre->restrictPoint(cell);
#endif
  size = _numQuadPts * _numBasis * _spaceDim;
  assert(size == _basisDerivPre->getFiberDimension(cell));
  for(int i=0; i < size; ++i)
    _basisDeriv[i] = values[i];
} // retrieveGeometry

// End of file 
