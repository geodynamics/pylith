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
typedef pylith::topology::SubMesh::SieveMesh SieveSubMesh;
typedef pylith::topology::SubMesh::RealSection SubRealSection;
typedef pylith::topology::Mesh::RealSection RealSection;
typedef pylith::topology::Mesh::RestrictVisitor RestrictVisitor;
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
  deallocate();
} // destructor

// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void 
pylith::faults::FaultCohesiveDyn::deallocate(void)
{ // deallocate
  FaultCohesive::deallocate();
} // deallocate
  
// ----------------------------------------------------------------------
// Initialize fault. Determine orientation and setup boundary
void
pylith::faults::FaultCohesiveDyn::initialize(const topology::Mesh& mesh,
					     const double upDir[3],
					     const double normalDir[3])
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
  assert(0 != _quadrature);
  assert(0 != _faultMesh); // changed to fault mesh
  assert(0 != _fields); // change to fields

  // Get cell geometry information that doesn't depend on cell
  const int numQuadPts = _quadrature->numQuadPts();
  const double_array& quadWts = _quadrature->quadWts();
  assert(quadWts.size() == numQuadPts);
  const int numBasis = _quadrature->numBasis();
  const int spaceDim = _quadrature->spaceDim();
  const int numCornersCohesive = 2*numBasis;

  // Allocate vectors for cell values.
  double_array tractionsCell(numQuadPts*spaceDim);
  double_array residualCell(numCornersCohesive*spaceDim);

  // GET COHESIVE CELLS (see FaultCohesiveKin)
  const ALE::Obj<SieveMesh>& sieveMesh = residual.mesh().sieveMesh();
  assert(!sieveMesh.isNull());
  const ALE::Obj<SieveMesh::label_sequence>& cellsCohesive = 
    sieveMesh->getLabelStratum("material-id", id());
  assert(!cellsCohesive.isNull());
  // Add setting cellsBegin

  // Get cell information
  const ALE::Obj<SieveSubMesh>& faultSieveMesh = _faultMesh->sieveMesh(); // updated (subSieveMesh -> faultSieveMesh)
  const SieveSubMesh::label_sequence::iterator cellsEnd = cellsCohesive->end();
  assert(!faultSieveMesh.isNull());

  // GET PARAMETERS FOR FAULT CONSTITUTIVE MODEL (later)
  // TEMPORARY hardwire to (1,0,0)
  for (int i=0; i < numQuadPts; ++i)
    tractionsCell[i*spaceDim] = 1.0;

  // Get sections (see FaultCohesiveKin)
  const ALE::Obj<RealSection>& residualSection = residual.section();
  topology::Mesh::UpdateAddVisitor residualVisitor(*residualSection,
						      &residualCell[0]);

#if !defined(PRECOMPUTE_GEOMETRY)
  double_array coordinatesCell(numBasis*spaceDim);
  const ALE::Obj<RealSection>& coordinates = 
    faultSieveMesh->getRealSection("coordinates");
  RestrictVisitor coordsVisitor(*coordinates,
				coordinatesCell.size(), &coordinatesCell[0]);
#endif

  // CHANGE TO LOOP OVER COHESIVE CELLS (See FaultCohesiveKin)
 
  // Loop over faces and integrate contribution from each face
  for (SieveSubMesh::label_sequence::iterator c_iter=cellsCohesive->begin();
       c_iter != cellsEnd; // cellsEnd
       ++c_iter) {
    const SieveMesh::point_type c_fault = _cohesiveToFault[*c_iter];
    residualCell = 0.0;

#if defined(PRECOMPUTE_GEOMETRY)
    _quadrature->retrieveGeometry(*c_iter);
#else
    coordsVisitor.clear();
    faultSieveMesh->restrictClosure(c_fault, coordsVisitor);
    _quadrature->computeGeometry(coordinatesCell, c_fault);
#endif

    // USE FAULT CONSTITUTIVE MODEL TO GET TRACTION (later)
    
    // Get cell geometry information that depends on cell
    const double_array& basis = _quadrature->basis();
    const double_array& jacobianDet = _quadrature->jacobianDet();

    // Compute action for traction bc terms
    for (int iQuad=0; iQuad < numQuadPts; ++iQuad) {
      const double wt = quadWts[iQuad] * jacobianDet[iQuad];
      for (int iBasis=0; iBasis < numBasis; ++iBasis) {
        const double valI = wt*basis[iQuad*numBasis+iBasis];
        for (int jBasis=0; jBasis < numBasis; ++jBasis) {
          const double valIJ = valI * basis[iQuad*numBasis+jBasis];
          for (int iDim=0; iDim < spaceDim; ++iDim)
            residualCell[iBasis*spaceDim+iDim] += 
	      tractionsCell[iQuad*spaceDim+iDim] * valIJ;
        } // for
      } // for
    } // for

    // Assemble cell contribution into field
    residualVisitor.clear();
    sieveMesh->updateAdd(*c_iter, residualVisitor);

    PetscLogFlops(numQuadPts*(1+numBasis*(1+numBasis*(1+2*spaceDim))));
  } // for
} // integrateResidual

// ----------------------------------------------------------------------
// Compute Jacobian matrix (A) associated with operator.
void
pylith::faults::FaultCohesiveDyn::integrateJacobian(
				   topology::Jacobian* jacobian,
				   const double t,
				   topology::SolutionFields* const fields)
{ // integrateJacobian
  // See Neumann
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
				       const topology::SolutionFields* fields)
{ // vertexField
  throw std::logic_error("FaultCohesiveDyn::vertexField() not implemented.");
} // vertexField

// ----------------------------------------------------------------------
// Get cell field associated with integrator.
const pylith::topology::Field<pylith::topology::SubMesh>&
pylith::faults::FaultCohesiveDyn::cellField(
				      const char* name,
				      const topology::SolutionFields* fields)
{ // cellField
  throw std::logic_error("FaultCohesiveDyn::cellField() not implemented.");
} // cellField

// ----------------------------------------------------------------------
// Calculate orientation at fault vertices.
void
pylith::faults::FaultCohesiveDyn::_calcOrientation(const double upDir[3],
						   const double normalDir[3])
{ // _calcOrientation
  throw std::logic_error("FaultCohesiveDyn::_calcOrientation() not implemented.");
} // _calcOrientation

// End of file 
