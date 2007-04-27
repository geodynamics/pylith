// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
//                           Brad T. Aagaard
//                        U.S. Geological Survey
//
// {LicenseText}
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "FaultCohesiveKin.hh" // implementation of object methods

#include "EqKinSrc.hh" // USES EqKinSrc

#include "pylith/feassemble/Quadrature.hh" // USES Quadrature
#include "pylith/utils/array.hh" // USES double_array

#include <assert.h> // USES assert()
#include <sstream> // USES std::ostringstream
#include <stdexcept> // USES std::runtime_error

// ----------------------------------------------------------------------
// Default constructor.
pylith::faults::FaultCohesiveKin::FaultCohesiveKin(void) :
  _eqsrc(0)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::faults::FaultCohesiveKin::~FaultCohesiveKin(void)
{ // destructor
  delete _eqsrc; _eqsrc = 0;
} // destructor

// ----------------------------------------------------------------------
// Copy constructor.
pylith::faults::FaultCohesiveKin::FaultCohesiveKin(const FaultCohesiveKin& f) :
  FaultCohesive(f),
  _eqsrc(0)
{ // copy constructor
  if (0 != f._eqsrc)
    _eqsrc = f._eqsrc->clone();
} // copy constructor

// ----------------------------------------------------------------------
// Set kinematic earthquake source.
void
pylith::faults::FaultCohesiveKin::eqsrc(EqKinSrc* src)
{ // eqsrc
  delete _eqsrc; _eqsrc = (0 != src) ? src->clone() : 0;
} // eqsrc

// ----------------------------------------------------------------------
// Initialize fault. Determine orientation and setup boundary
void
pylith::faults::FaultCohesiveKin::initialize(const ALE::Obj<ALE::Mesh>& mesh,
					     const double_array& upDir)
{ // initialize
  assert(0 != _quadrature);
  assert(0 != _faultMesh);
  assert(!_faultMesh->isNull());
  
  if (3 != upDir.size())
    throw std::runtime_error("Up direction for fault orientation must be "
			     "a vector with 3 components.");

  // Allocate section for orientation at all fault vertices
  ALE::Obj<real_section_type> orientation = 
    new real_section_type((*_faultMesh)->comm(), (*_faultMesh)->debug());
  assert(!orientation.isNull());
  const int orientationSize = _orientationSize();
  orientation->setFiberDimension((*_faultMesh)->depthStratum(0), 
				 orientationSize);
  (*_faultMesh)->allocate(orientation);
  
  // Get section containing coordinates of vertices
  const ALE::Obj<real_section_type>& coordinates = 
    mesh->getRealSection("coordinates");
  assert(!coordinates.isNull());

  // Loop over cells
  const ALE::Obj<Mesh::label_sequence>& cells = 
    (*_faultMesh)->heightStratum(0);
  const Mesh::label_sequence::iterator cBegin = cells->begin();
  const Mesh::label_sequence::iterator cEnd = cells->end();
  double_array cellOrientation(_quadrature->numBasis()*orientationSize);
  for (Mesh::label_sequence::iterator c_iter=cBegin;
       c_iter != cEnd;
       ++c_iter) {
    // Compute cell geometry at vertices
    _quadrature->computeGeometryVert(*_faultMesh, coordinates, *c_iter);

    const double_array& jacobian = _quadrature->jacobianVert();
    const double_array& jacobianDet = _quadrature->jacobianDetVert();

    // Compute weighted orientation of face at vertices (using geometry info)
    // NEED TO CALL APPROPRIATE ORIENT FOR CELL DIM
    //_orient1D(&cellOrientation, jacobian, jacobianDet, upDir, numVertices);

    // Update weighted orientations
    // Loop over vertices in cell
    //   Update section with orientation for each vertex
  } // for

  // Assemble orientation information

  // Loop over vertices
  //   Make orientation information unit magnitude

  // Create list of constraint vertices

  // Only store orientation information at constraint vertices
} // initialize

// ----------------------------------------------------------------------
// Integrate contribution of cohesive cells to residual term.
void
pylith::faults::FaultCohesiveKin::integrateResidual(
				const ALE::Obj<real_section_type>& residual,
				const ALE::Obj<real_section_type>& disp,
				const ALE::Obj<Mesh>& mesh)
{ // integrateResidual

  // Subtract constraint forces (which are in disp at the constraint
  // DOF) to residual; contributions are at DOF of normal vertices (i and j)

} // integrateResidual

// ----------------------------------------------------------------------
// Compute Jacobian matrix (A) associated with operator.
void
pylith::faults::FaultCohesiveKin::integrateJacobian(
				    PetscMat* mat,
				    const ALE::Obj<real_section_type>& dispT,
				    const ALE::Obj<Mesh>& mesh)
{ // integrateJacobian

  // Add constraint information to Jacobian matrix; these are the
  // direction cosines. Entries are associated with vertices ik, jk,
  // ki, and kj.

} // integrateJacobian
  
// ----------------------------------------------------------------------
// Set field.
void
pylith::faults::FaultCohesiveKin::setField(
				     const double t,
				     const ALE::Obj<real_section_type>& disp,
				     const ALE::Obj<Mesh>& mesh)
{ // setField
  typedef std::vector<Mesh::point_type>::const_iterator vert_iterator;

  assert(0 != _eqsrc);

  const ALE::Obj<real_section_type>& slip = _eqsrc->slip(t, _constraintVert);
  assert(!slip.isNull());
  const vert_iterator vBegin = _constraintVert.begin();
  const vert_iterator vEnd = _constraintVert.end();
  for (vert_iterator v_iter=vBegin; v_iter != vEnd; ++v_iter)
    disp->updatePoint(*v_iter, slip->restrictPoint(*v_iter));
} // setField

// ----------------------------------------------------------------------
// Get size (fiber dimension) of orientation information.
int 
pylith::faults::FaultCohesiveKin::_orientationSize(void) const
{ // _orientationSize
  assert(0 != _quadrature);

  int size = 0;
  switch (_quadrature->cellDim()) {
  case 1 :
    size = 1;
    break;
  case 2 :
    size = 4;
    break;
  case 3 :
    size = 9;
    break;
  default :
    assert(0);
  } // switch
  return size;
} // _orientationSize


// End of file 
