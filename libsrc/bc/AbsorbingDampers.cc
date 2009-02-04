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

#include "AbsorbingDampers.hh" // implementation of object methods

#include "pylith/feassemble/Quadrature.hh" // USES Quadrature
#include "pylith/feassemble/CellGeometry.hh" // USES CellGeometry
#include "pylith/topology/SolutionFields.hh" // USES SolutionFields

#include "spatialdata/spatialdb/SpatialDB.hh" // USES SpatialDB
#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include <cstring> // USES memcpy()
#include <cassert> // USES assert()
#include <stdexcept> // USES std::runtime_error
#include <sstream> // USES std::ostringstream

// ----------------------------------------------------------------------
// Default constructor.
pylith::bc::AbsorbingDampers::AbsorbingDampers(void) :
  _boundaryMesh(0),
  _dampingConsts(0)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::bc::AbsorbingDampers::~AbsorbingDampers(void)
{ // destructor
  delete _boundaryMesh; _boundaryMesh = 0;
  delete _dampingConsts; _dampingConsts = 0;
} // destructor

// ----------------------------------------------------------------------
// Initialize boundary condition. Determine orienation and compute traction
// vector at integration points.
void
pylith::bc::AbsorbingDampers::initialize(const topology::Mesh& mesh,
					 const double upDir[3])
{ // initialize
  assert(0 != _quadrature);
  assert(0 != _db);

  delete _boundaryMesh; _boundaryMesh = 0;
  delete _dampingConst; _dampingConsts = 0;

  _boundaryMesh = new topology::SubMesh(mesh, _label.c_str());
  assert(0 != _boundaryMesh);

  // check compatibility of quadrature and boundary mesh
  if (_quadrature->cellDim() != _boundaryMesh->dimension()) {
    std::ostringstream msg;
    msg << "Quadrature is incompatible with cells for absorbing boundary "
	<< "condition '" << _label << "'.\n"
	<< "Dimension of boundary mesh: " << _boundaryMesh->dimension()
	<< ", dimension of quadrature: " << _quadrature->cellDim()
	<< ".";
    throw std::runtime_error(msg.str());
  } // if
  const int numCorners = _quadrature->numBasis();

  // Get 'surface' cells (1 dimension lower than top-level cells)
  const ALE::Obj<SieveSubMesh>& submesh = _boundaryMesh->sieveMesh();
  assert(!submesh.isNull());
  const ALE::Obj<SubMesh::label_sequence>& cells = submesh->heightStratum(1);
  assert(!cells.isNull());
  const SieveSubMesh::label_sequence::iterator cellsEnd = cells->end();
  const int boundaryDepth = submesh->depth()-1; // depth of bndry cells
  for (SieveSubMesh::label_sequence::iterator c_iter=cells->begin();
       c_iter != cellsEnd;
       ++c_iter) {
    const int cellNumCorners = 
      submesh->getNumCellCorners(*c_iter, boundaryDepth);

    if (numCorners != cellNumCorners) {
      std::ostringstream msg;
      msg << "Quadrature is incompatible with cell for absorbing boundary "
	  << "condition '" << _label << "'.\n"
	  << "Cell " << *c_iter << " has " << cellNumCorners
	  << " vertices but quadrature reference cell has "
	  << numCorners << " vertices.";
      throw std::runtime_error(msg.str());
    } // if
  } // for

  // Get damping constants at each quadrature point and rotate to
  // global coordinate frame using orientation information
  const int cellDim = _quadrature->cellDim();
  const int numBasis = _quadrature->numBasis();
  const int numQuadPts = _quadrature->numQuadPts();
  const int spaceDim = cs->spaceDim();
  const int fiberDim = numQuadPts * spaceDim;
  _dampingConsts = new FieldSubMesh(*_boundaryMesh);
  assert(0 != _dampingConsts);
  _dampingConsts->newSection(cells, fiberDim);
  _dampingConsts->allocate();

  // Containers for orientation information
  const int orientationSize = spaceDim*spaceDim;
  const feassemble::CellGeometry& cellGeometry = _quadrature->refGeometry();
  const int jacobianSize = (cellDim > 0) ? spaceDim * cellDim : 1;
  double_array jacobian(jacobianSize);
  double jacobianDet = 0;
  double_array orientation(orientationSize);
  double_array cellVertices(numBasis*spaceDim);

  // open database with material property information
  _db->open();
  int numValues = 0;
  if (cellDim > 0) {
    numValues = 3;
    const char* valueNames[] = { "density", "vp", "vs" };
    _db->queryVals(valueNames, numValues);
  } else {
    numValues = 2;
    const char* valueNames[] = { "density", "vp" };
    _db->queryVals(valueNames, numValues);
  } // else

  // Container for data returned in query of database
  double_array queryData(numValues);
  double_array quadPtRef(cellDim);
  double_array quadPtsGlobal(numQuadPts*spaceDim);

  // Container for damping constants for current cell
  double_array dampingConstsLocal(fiberDim);
  double_array dampingConstsGlobal(fiberDim);

  const ALE::Obj<MeshRealSection>& coordinates =
    mesh->getRealSection("coordinates");
  assert(!coordinates.isNull());

  assert(0 != _normalizer);
  const double lengthScale = _normalizer->lengthScale();
  const double densityScale = _normalizer->densityScale();
  const double velocityScale = 
    _normalizer->lengthScale() / _normalizer->timeScale();

  const ALE::Obj<SubMeshRealSection>& constsSection = 
    _dampingConsts->section();
  assert(!constsSection.isNull());

  for(SieveSubMesh::label_sequence::iterator c_iter = cells->begin();
      c_iter != cellsEnd;
      ++c_iter) {
    _quadrature->computeGeometry(submesh, coordinates, *c_iter);
    const double_array& quadPtsNondim = _quadrature->quadPts();
    const double_array& quadPtsRef = _quadrature->quadPtsRef();
    quadPtsGlobal = quadPtsNondim;
    _normalizer->dimensionalize(&quadPtsGlobal[0], quadPtsGlobal.size(), 
				lengthScale);

    dampingConstsGlobal = 0.0;
    for(int iQuad = 0; iQuad < numQuadPts; ++iQuad) {
      // Compute damping constants in normal/tangential coordinates
      const int err = _db->query(&queryData[0], numValues, 
				 &quadPtsGlobal[iQuad*spaceDim], spaceDim, cs);
      if (err) {
	std::ostringstream msg;
	msg << "Could not find parameters for physical properties at \n"
	    << "(";
	for (int i=0; i < spaceDim; ++i)
	  msg << "  " << quadPtsGlobal[iQuad*spaceDim+i];
	msg << ") for absorbing boundary condition " << _label << "\n"
	    << "using spatial database " << _db->label() << ".";
	throw std::runtime_error(msg.str());
      } // if
      // Nondimensionalize damping constants
      const double densityN = 
	_normalizer->nondimensionalize(queryData[0], densityScale);
      const double vpN = 
	_normalizer->nondimensionalize(queryData[1], velocityScale);
      const double vsN = (3 == numValues) ?
	_normalizer->nondimensionalize(queryData[2], velocityScale) :
	0.0;
      
      const double constTangential = densityN * vsN;
      const double constNormal = densityN * vpN;
      const int numTangential = spaceDim-1;
      for (int iDim=0; iDim < numTangential; ++iDim)
	dampingConstsLocal[iDim] = constTangential;
      dampingConstsLocal[spaceDim-1] = constNormal;

      // Compute normal/tangential orientation
      submesh->restrictClosure(coordinates, *c_iter, 
			       &cellVertices[0], cellVertices.size());
      memcpy(&quadPtRef[0], &quadPtsRef[iQuad*cellDim], cellDim*sizeof(double));
      cellGeometry.jacobian(&jacobian, &jacobianDet, cellVertices, quadPtRef);
      cellGeometry.orientation(&orientation, jacobian, jacobianDet, 
			       upDir);
      orientation /= jacobianDet;

      for (int iDim=0; iDim < spaceDim; ++iDim) {
	for (int jDim=0; jDim < spaceDim; ++jDim)
	  dampingConstsGlobal[iQuad*spaceDim+iDim] += 
	    dampingConstsLocal[jDim]*orientation[jDim*spaceDim+iDim];
	// Ensure damping constants are positive
	dampingConstsGlobal[iQuad*spaceDim+iDim] = 
	  fabs(dampingConstsGlobal[iQuad*spaceDim+iDim]);
      } // for
    } // for
    constsSection->updatePoint(*c_iter, &dampingConstsGlobal[0]);
  } // for

  _db->close();
} // initialize

// ----------------------------------------------------------------------
// Integrate contributions to residual term (r) for operator.
void
pylith::bc::AbsorbingDampers::integrateResidual(
				   const topology::Field& residual,
				 const double t,
				 topology::SolutionFields* const fields)
{ // integrateResidual
  assert(0 != _quadrature);
  assert(0 != _boundaryMesh);
  assert(0 != fields);

  PetscErrorCode err = 0;

  // Get cell information
  const ALE::Obj<SieveSubMesh>& submesh = _boundaryMesh->sieveMesh();
  assert(!submesh.isNull());
  const ALE::Obj<SubMesh::label_sequence>& cells = submesh->heightStratum(1);
  assert(!cells.isNull());
  const SubMesh::label_sequence::iterator cellsEnd = cells->end();

  // Get sections
  const ALE::Obj<MeshRealSection>& coordinates = 
    mesh->getRealSection("coordinates");
  assert(!coordinates.isNull());
  const ALE::Obj<MeshRealSection>& dispTpdt = fields->getField("disp t+dt");
  const ALE::Obj<MeshRealSection>& dispTmdt = fields->getField("disp t-dt");
  assert(!dispTpdt.isNull());
  assert(!dispTmdt.isNull());

  // Get parameters used in integration.
  const double dt = _dt;
  assert(dt > 0);

  // Get cell geometry information that doesn't depend on cell
  const int numQuadPts = _quadrature->numQuadPts();
  const double_array& quadWts = _quadrature->quadWts();
  assert(quadWts.size() == numQuadPts);
  const int numBasis = _quadrature->numBasis();
  const int spaceDim = _quadrature->spaceDim();
  const int cellDim = _quadrature->cellDim();

  // Allocate vectors for cell values.
  _initCellVector();
  const int cellVecSize = numBasis*spaceDim;
  double_array dispTpdtCell(cellVecSize);
  double_array dispTmdtCell(cellVecSize);

  for (SieveSubMesh::label_sequence::iterator c_iter=cells->begin();
       c_iter != cellsEnd;
       ++c_iter) {
    // Compute geometry information for current cell
    _quadrature->computeGeometry(submesh, coordinates, *c_iter);

    // Reset element vector to zero
    _resetCellVector();

    // Restrict input fields to cell
    submesh->restrictClosure(dispTpdt, *c_iter, 
			     &dispTpdtCell[0], cellVecSize);
    submesh->restrictClosure(dispTmdt, *c_iter, 
			     &dispTmdtCell[0], cellVecSize);
    assert(numQuadPts*spaceDim == _dampingConsts->getFiberDimension(*c_iter));
    const real_section_type::value_type* dampingConstsCell = 
      _dampingConsts->restrictPoint(*c_iter);

    // Get cell geometry information that depends on cell
    const double_array& basis = _quadrature->basis();
    const double_array& jacobianDet = _quadrature->jacobianDet();

    // Compute action for absorbing bc terms
    for (int iQuad=0; iQuad < numQuadPts; ++iQuad) {
      const double wt = 
	quadWts[iQuad] * jacobianDet[iQuad] / (2.0 * dt);

      for (int iBasis=0; iBasis < numBasis; ++iBasis) {
        const double valI = wt*basis[iQuad*numBasis+iBasis];
        for (int jBasis=0; jBasis < numBasis; ++jBasis) {
          const double valIJ = valI * basis[iQuad*numBasis+jBasis];
          for (int iDim=0; iDim < spaceDim; ++iDim)
            _cellVector[iBasis*spaceDim+iDim] += 
	      dampingConstsCell[iQuad*spaceDim+iDim] *
	      valIJ * (dispTmdtCell[jBasis*spaceDim+iDim]);
        } // for
      } // for

    } // for
    PetscLogFlops(numQuadPts*(3+numBasis*(1+numBasis*(5*spaceDim))));

    // Assemble cell contribution into field
    submesh->updateAdd(residual, *c_iter, _cellVector);
  } // for
} // integrateResidual

// ----------------------------------------------------------------------
// Integrate contributions to Jacobian matrix (A) associated with
void
pylith::bc::AbsorbingDampers::integrateJacobian(
					 PetscMat* jacobian,
					 const double t,
					 topology::SolutionFields* const fields)
{ // integrateJacobian
  assert(0 != _quadrature);
  assert(0 != _boundaryMesh);
  assert(0 != jacobian);
  assert(0 != fields);

  typedef ALE::ISieveVisitor::IndicesVisitor<MeshRealSection,SieveMesh::order_type,PetscInt> visitor_type;
  PetscErrorCode err = 0;

  // Get cell information
  const ALE::Obj<SieveSubMesh>& submesh = _boundaryMesh->sieveMesh();
  assert(!submesh.isNull());
  const ALE::Obj<SieveSubMesh::label_sequence>& cells = 
    submesh->heightStratum(1);
  assert(!cells.isNull());
  const SieveSubMesh::label_sequence::iterator cellsEnd = cells->end();

  // Get sections
  const ALE::Obj<MeshRealSection>& coordinates = 
    mesh->getRealSection("coordinates");
  assert(!coordinates.isNull());
  const ALE::Obj<MeshRealSection>& dispT = fields->getField("disp t");
  assert(!dispT.isNull());

  // Get parameters used in integration.
  const double dt = _dt;
  assert(dt > 0);

  // Get cell geometry information that doesn't depend on cell
  const int numQuadPts = _quadrature->numQuadPts();
  const double_array& quadWts = _quadrature->quadWts();
  assert(quadWts.size() == numQuadPts);
  const int numBasis = _quadrature->numBasis();
  const int spaceDim = _quadrature->spaceDim();
  const int cellDim = _quadrature->cellDim();

  // Allocate vectors for cell values.
  _initCellVector();
  const int cellVecSize = numBasis*spaceDim;

  // Allocate vector for cell values (if necessary)
  _initCellMatrix();

  const ALE::Obj<MeshRealSection>& solution = fields->getSolution();
  assert(!solution.isNull());  
  const ALE::Obj<SieveMesh::order_type>& globalOrder = 
    mesh->getFactory()->getGlobalOrder(mesh, "default", solution);
  assert(!globalOrder.isNull());
  visitor_type iV(*solution, *globalOrder,
		  (int) pow(mesh->getSieve()->getMaxConeSize(),
			    mesh->depth())*spaceDim);

  for (SieveSubMesh::label_sequence::iterator c_iter=cells->begin();
       c_iter != cellsEnd;
       ++c_iter) {
    // Compute geometry information for current cell
    _quadrature->computeGeometry(_boundaryMesh, coordinates, *c_iter);

    // Reset element vector to zero
    _resetCellMatrix();

    // Get cell geometry information that depends on cell
    const double_array& basis = _quadrature->basis();
    const double_array& jacobianDet = _quadrature->jacobianDet();

    // Restrict input fields to cell
    assert(numQuadPts*spaceDim == _dampingConsts->getFiberDimension(*c_iter));
    const MeshRealSection::value_type* dampingConstsCell = 
      _dampingConsts->restrictPoint(*c_iter);

    // Compute Jacobian for absorbing bc terms
    for (int iQuad=0; iQuad < numQuadPts; ++iQuad) {
      const double wt = 
	quadWts[iQuad] * jacobianDet[iQuad] / (2.0 * dt);
      for (int iBasis=0, iQ=iQuad*numBasis; iBasis < numBasis; ++iBasis) {
        const double valI = wt*basis[iQ+iBasis];
        for (int jBasis=0; jBasis < numBasis; ++jBasis) {
          const double valIJ = valI * basis[iQ+jBasis];
          for (int iDim=0; iDim < spaceDim; ++iDim) {
            const int iBlock = (iBasis*spaceDim + iDim) * (numBasis*spaceDim);
            const int jBlock = (jBasis*spaceDim + iDim);
            _cellMatrix[iBlock+jBlock] += 
	      valIJ * dampingConstsCell[iQuad*spaceDim+iDim];
          } // for
        } // for
      } // for
    } // for
    PetscLogFlops(numQuadPts*(3+numBasis*(1+numBasis*(1+2*spaceDim))));
    
    // Assemble cell contribution into PETSc Matrix
    err = updateOperator(*jacobian, submesh->getSieve(), 
			 iV, *c_iter, _cellMatrix, ADD_VALUES);
    if (err)
      throw std::runtime_error("Update to PETSc Mat failed.");
    iV.clear();
  } // for

  _needNewJacobian = false;
} // integrateJacobian

// ----------------------------------------------------------------------
// Verify configuration is acceptable.
void
pylith::bc::AbsorbingDampers::verifyConfiguration(const topology::Mesh& mesh) const
{ // verifyConfiguration
  BoundaryCondition::verifyConfiguration(mesh);
} // verifyConfiguration


// End of file 
