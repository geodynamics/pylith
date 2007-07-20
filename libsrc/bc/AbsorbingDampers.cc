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
#include "pylith/topology/FieldsManager.hh" // USES FieldsManager

#include "spatialdata/spatialdb/SpatialDB.hh" // USES SpatialDB
#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys

#include <Selection.hh> // USES submesh algorithms

#include <assert.h> // USES assert()
#include <stdexcept> // USES std::runtime_error
#include <sstream> // USES std::ostringstream

// ----------------------------------------------------------------------
// Default constructor.
pylith::bc::AbsorbingDampers::AbsorbingDampers(void)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::bc::AbsorbingDampers::~AbsorbingDampers(void)
{ // destructor
} // destructor

// ----------------------------------------------------------------------
// Initialize boundary condition. Determine orienation and compute traction
// vector at integration points.
void
pylith::bc::AbsorbingDampers::initialize(const ALE::Obj<ALE::Mesh>& mesh,
				const spatialdata::geocoords::CoordSys* cs,
				const double_array& upDir)
{ // initialize
  assert(0 != _quadrature);
  assert(0 != _db);
  assert(!mesh.isNull());
  assert(0 != cs);

  if (3 != upDir.size())
    throw std::runtime_error("Up direction for boundary orientation must be "
			     "a vector with 3 components.");

  // Extract submesh associated with boundary
  _boundaryMesh = 
    ALE::Selection<ALE::Mesh>::submesh(mesh, mesh->getIntSection(_label));
  if (_boundaryMesh.isNull()) {
    std::ostringstream msg;
    msg << "Could not construct boundary mesh for absorbing boundary "
	<< "condition '" << _label << "'.";
    throw std::runtime_error(msg.str());
  } // if

  // check compatibility of quadrature and boundary mesh
  if (_quadrature->cellDim() != _boundaryMesh->getDimension()) {
    std::ostringstream msg;
    msg << "Quadrature is incompatible with cells for absorbing boundary "
	<< "condition '" << _label << "'.\n"
	<< "Dimension of boundary mesh: " << _boundaryMesh->getDimension()
	<< ", dimension of quadrature: " << _quadrature->cellDim()
	<< ".";
    throw std::runtime_error(msg.str());
  } // if
  const int numCorners = _quadrature->numBasis();

  // Get 'surface' cells (1 dimension lower than top-level cells)
  const ALE::Obj<ALE::Mesh::label_sequence>& cells = 
    _boundaryMesh->heightStratum(1);
  assert(!cells.isNull());
  const Mesh::label_sequence::iterator cellsBegin = cells->begin();
  const Mesh::label_sequence::iterator cellsEnd = cells->end();
  const ALE::Obj<sieve_type>& sieve = mesh->getSieve();
  assert(!sieve.isNull());
  for (Mesh::label_sequence::iterator c_iter=cellsBegin;
       c_iter != cellsEnd;
       ++c_iter) {
    const int cellNumCorners = sieve->nCone(*c_iter, mesh->depth())->size();
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
  const int fiberDim = spaceDim * numQuadPts;
  _dampingConsts = new real_section_type(mesh->comm(), mesh->debug());
  assert(!_dampingConsts.isNull());
  _dampingConsts->setFiberDimension(cells, fiberDim);
  mesh->allocate(_dampingConsts);

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
  if (cellDim > 1) {
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
  double_array quadPtRef(spaceDim);

  // Container for damping constants for current cell
  double_array dampingConstsLocal(fiberDim);
  double_array dampingConstsGlobal(fiberDim);

  const ALE::Obj<real_section_type>& coordinates =
    mesh->getRealSection("coordinates");

  for(Mesh::label_sequence::iterator c_iter = cells->begin();
      c_iter != cells->end();
      ++c_iter) {
    _quadrature->computeGeometry(_boundaryMesh, coordinates, *c_iter);
    const double_array& quadPts = _quadrature->quadPts();
    const double_array& quadPtsRef = _quadrature->quadPtsRef();

    for(int iQuad = 0, index=0; iQuad < numQuadPts; ++iQuad, index+=spaceDim) {
      // Compute damping constants in normal/tangential coordinates
      const int err = _db->query(&queryData[0], numValues, 
				 &quadPts[index], spaceDim, cs);
      if (err) {
	std::ostringstream msg;
	msg << "Could not find parameters for physical properties at \n"
	    << "(";
	for (int i=0; i < spaceDim; ++i)
	  msg << "  " << quadPts[index+spaceDim];
	msg << ") for absorbing boundary condition " << _label << "\n"
	    << "using spatial database " << _db->label() << ".";
	throw std::runtime_error(msg.str());
      } // if
      dampingConstsLocal[index  ] = queryData[0]*queryData[1];
      for (int iDim=1; iDim < spaceDim; ++iDim)
	dampingConstsLocal[index+iDim] = queryData[0]*queryData[2];

      // Compute normal/tangential orientation
      mesh->restrict(coordinates, *c_iter, 
		     &cellVertices[0], cellVertices.size());
      memcpy(&quadPtRef[0], &quadPtsRef[index], spaceDim*sizeof(double));
      cellGeometry.jacobian(&jacobian, &jacobianDet, cellVertices, quadPtRef);
      cellGeometry.orientation(&orientation, jacobian, jacobianDet, 
			       upDir);
      dampingConstsGlobal = 0.0;
      for (int iDim=0; iDim < spaceDim; ++iDim)
	for (int jDim=0; jDim < spaceDim; ++jDim)
	  dampingConstsGlobal[iDim] += 
	    dampingConstsLocal[jDim]*orientation[iDim*spaceDim+jDim];
    } // for
    mesh->update(_dampingConsts, *c_iter, &dampingConstsGlobal[0]);
  } // for

  _db->close();
} // initialize

// ----------------------------------------------------------------------
// Integrate contributions to residual term (r) for operator.
void
pylith::bc::AbsorbingDampers::integrateResidual(
				 const ALE::Obj<real_section_type>& residual,
				 const double t,
				 topology::FieldsManager* const fields,
				 const ALE::Obj<Mesh>& mesh)
{ // integrateResidual
  assert(0 != _quadrature);
  assert(!_boundaryMesh.isNull());
  assert(!residual.isNull());
  assert(0 != fields);
  assert(!mesh.isNull());

  PetscErrorCode err = 0;

  // Get cell information
  const ALE::Obj<ALE::Mesh::label_sequence>& cells = 
    _boundaryMesh->heightStratum(1);
  assert(!cells.isNull());
  const Mesh::label_sequence::iterator cellsEnd = cells->end();

  // Get sections
  const ALE::Obj<real_section_type>& coordinates = 
    mesh->getRealSection("coordinates");
  assert(!coordinates.isNull());
  const ALE::Obj<real_section_type>& dispTpdt = fields->getHistoryItem(0);
  const ALE::Obj<real_section_type>& dispTmdt = fields->getHistoryItem(2);
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
  double_array dampingConstsCell(numQuadPts*spaceDim);

  for (Mesh::label_sequence::iterator c_iter=cells->begin();
       c_iter != cellsEnd;
       ++c_iter) {
    // Compute geometry information for current cell
    _quadrature->computeGeometry(_boundaryMesh, coordinates, *c_iter);

    // Reset element vector to zero
    _resetCellVector();

    // Restrict input fields to cell
    _boundaryMesh->restrict(dispTpdt, *c_iter, &dispTpdtCell[0], cellVecSize);
    _boundaryMesh->restrict(dispTmdt, *c_iter, &dispTmdtCell[0], cellVecSize);
    _boundaryMesh->restrict(_dampingConsts, *c_iter, 
			    &dampingConstsCell[0], dampingConstsCell.size());

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
	      valIJ * (- dispTpdtCell[jBasis*spaceDim+iDim] 
		       + dispTmdtCell[jBasis*spaceDim+iDim]);
        } // for
      } // for
    } // for
    err = PetscLogFlops(numQuadPts*(3+numBasis*(1+numBasis*(5*spaceDim))));
    if (err)
      throw std::runtime_error("Logging PETSc flops failed.");

    // Assemble cell contribution into field
    mesh->updateAdd(residual, *c_iter, _cellVector);
  } // for
} // integrateResidual

// ----------------------------------------------------------------------
// Integrate contributions to Jacobian matrix (A) associated with
void
pylith::bc::AbsorbingDampers::integrateJacobian(
					 PetscMat* jacobian,
					 const double t,
					 topology::FieldsManager* const fields,
					 const ALE::Obj<Mesh>& mesh)
{ // integrateJacobian
  assert(0 != _quadrature);
  assert(!_boundaryMesh.isNull());
  assert(0 != jacobian);
  assert(0 != fields);
  assert(!mesh.isNull());

  PetscErrorCode err = 0;

  // Get cell information
  const ALE::Obj<ALE::Mesh::label_sequence>& cells = 
    _boundaryMesh->heightStratum(1);
  assert(!cells.isNull());
  const Mesh::label_sequence::iterator cellsEnd = cells->end();

  // Get sections
  const ALE::Obj<real_section_type>& coordinates = 
    mesh->getRealSection("coordinates");
  assert(!coordinates.isNull());
  const ALE::Obj<real_section_type>& dispT = fields->getHistoryItem(1);
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
  double_array dampingConstsCell(numQuadPts*spaceDim);

  // Allocate vector for cell values (if necessary)
  _initCellMatrix();

  for (Mesh::label_sequence::iterator c_iter=cells->begin();
       c_iter != cellsEnd;
       ++c_iter) {
    // Compute geometry information for current cell
    _quadrature->computeGeometry(_boundaryMesh, coordinates, *c_iter);

    // Reset element vector to zero
    _resetCellMatrix();

    // Get cell geometry information that depends on cell
    const double_array& basis = _quadrature->basis();
    const double_array& jacobianDet = _quadrature->jacobianDet();

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
    PetscErrorCode err = 
      PetscLogFlops(numQuadPts*(3+numBasis*(1+numBasis*(1+2*spaceDim))));
    if (err)
      throw std::runtime_error("Logging PETSc flops failed.");
    
    // Assemble cell contribution into PETSc Matrix
    const ALE::Obj<Mesh::order_type>& globalOrder = 
      mesh->getFactory()->getGlobalOrder(mesh, "default", dispT);
    assert(!globalOrder.isNull());

    err = updateOperator(*jacobian, mesh, dispT, globalOrder,
			 *c_iter, _cellMatrix, ADD_VALUES);
    if (err)
      throw std::runtime_error("Update to PETSc Mat failed.");
  } // for

  _needNewJacobian = false;
} // integrateJacobian

// ----------------------------------------------------------------------
// Verify configuration is acceptable.
void
pylith::bc::AbsorbingDampers::verifyConfiguration(const ALE::Obj<Mesh>& mesh)
{ // verifyConfiguration
  throw std::logic_error("Not implemented.");
} // verifyConfiguration


// End of file 
