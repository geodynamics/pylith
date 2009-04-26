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

#include "FaultCohesiveDyn.hh" // implementation of object methods

#include <cassert> // USES assert()
#include <sstream> // USES std::ostringstream
#include <stdexcept> // USES std::runtime_error

// ----------------------------------------------------------------------
typedef pylith::topology::Mesh::SieveMesh SieveMesh;

// ----------------------------------------------------------------------
// Default constructor.
pylith::faults::FaultCohesiveDyn::FaultCohesiveDyn(void)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::faults::FaultCohesiveDyn::~FaultCohesiveDyn(void)
{ // destructor
} // destructor

// ----------------------------------------------------------------------
// Initialize fault. Determine orientation and setup boundary
void
pylith::faults::FaultCohesiveDyn::initialize(const topology::Mesh& mesh,
					     const double upDir[3],
					     const double normalDir[3],
					     spatialdata::spatialdb::SpatialDB* matDB)
{ // initialize
  throw std::logic_error("FaultCohesiveDyn::initialize() not implemented.");
} // initialize

// ----------------------------------------------------------------------
// Integrate contribution of cohesive cells to residual term.
void
pylith::faults::FaultCohesiveDyn::integrateResidual(
			   const topology::Field<topology::Mesh>& residual,
			   const double t,
			   topology::SolutionFields* const fields)
{ // integrateResidual
  throw std::logic_error("FaultCohesiveDyn::integrateResidual() not implemented.");
} // integrateResidual

// ----------------------------------------------------------------------
// Compute Jacobian matrix (A) associated with operator.
void
pylith::faults::FaultCohesiveDyn::integrateJacobian(
				   topology::Jacobian* jacobian,
				   const double t,
				   topology::SolutionFields* const fields)
{ // integrateJacobian
  throw std::logic_error("FaultCohesiveDyn::integrateJacobian() not implemented.");
} // integrateJacobian
  
// ----------------------------------------------------------------------
// Verify configuration is acceptable.
void
pylith::faults::FaultCohesiveDyn::verifyConfiguration(
					    const topology::Mesh& mesh) const
{ // verifyConfiguration
  assert(0 != _quadrature);

  const ALE::Obj<SieveMesh>& sieveMesh = mesh.sieveMesh();
  assert(!sieveMesh.isNull());

  if (!sieveMesh->hasIntSection(label())) {
    std::ostringstream msg;
    msg << "Mesh missing group of vertices '" << label()
	<< " for boundary condition.";
    throw std::runtime_error(msg.str());
  } // if  

  // check compatibility of mesh and quadrature scheme
  const int dimension = mesh.dimension()-1;
  if (_quadrature->cellDim() != dimension) {
    std::ostringstream msg;
    msg << "Dimension of reference cell in quadrature scheme (" 
	<< _quadrature->cellDim() 
	<< ") does not match dimension of cells in mesh (" 
	<< dimension << ") for fault '" << label()
	<< "'.";
    throw std::runtime_error(msg.str());
  } // if
  const int numCorners = _quadrature->numBasis();
  const ALE::Obj<SieveMesh::label_sequence>& cells = 
    sieveMesh->getLabelStratum("material-id", id());
  assert(!cells.isNull());
  const SieveMesh::label_sequence::iterator cellsBegin = cells->begin();
  const SieveMesh::label_sequence::iterator cellsEnd = cells->end();
  for (SieveMesh::label_sequence::iterator c_iter=cellsBegin;
       c_iter != cellsEnd;
       ++c_iter) {
    const int cellNumCorners = sieveMesh->getNumCellCorners(*c_iter);
    if (3*numCorners != cellNumCorners) {
      std::ostringstream msg;
      msg << "Number of vertices in reference cell (" << numCorners 
	  << ") is not compatible with number of vertices (" << cellNumCorners
	  << ") in cohesive cell " << *c_iter << " for fault '"
	  << label() << "'.";
      throw std::runtime_error(msg.str());
    } // if
  } // for
} // verifyConfiguration

// ----------------------------------------------------------------------
// Get vertex field associated with integrator.
const pylith::topology::Field<pylith::topology::SubMesh>&
pylith::faults::FaultCohesiveDyn::vertexField(
				       const char* name,
				       const topology::SolutionFields& fields)
{ // vertexField
  throw std::logic_error("FaultCohesiveDyn::vertexField() not implemented.");
} // vertexField

// ----------------------------------------------------------------------
// Get cell field associated with integrator.
const pylith::topology::Field<pylith::topology::SubMesh>&
pylith::faults::FaultCohesiveDyn::cellField(
				      const char* name,
				      const topology::SolutionFields& fields)
{ // cellField
  throw std::logic_error("FaultCohesiveDyn::vertexField() not implemented.");
} // cellField


// End of file 
