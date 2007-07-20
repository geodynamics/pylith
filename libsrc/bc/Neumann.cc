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

#include "Neumann.hh" // implementation of object methods

#include "pylith/feassemble/Quadrature.hh" // USES Quadrature
#include "pylith/feassemble/CellGeometry.hh" // USES CellGeometry

#include "spatialdata/spatialdb/SpatialDB.hh" // USES SpatialDB
#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys

#include <Selection.hh> // USES submesh algorithms

#include <assert.h> // USES assert()
#include <stdexcept> // USES std::runtime_error
#include <sstream> // USES std::ostringstream

// ----------------------------------------------------------------------
// Default constructor.
pylith::bc::Neumann::Neumann(void)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::bc::Neumann::~Neumann(void)
{ // destructor
} // destructor

// ----------------------------------------------------------------------
// Initialize boundary condition. Determine orienation and compute traction
// vector at integration points.
void
pylith::bc::Neumann::initialize(const ALE::Obj<ALE::Mesh>& mesh,
				const spatialdata::geocoords::CoordSys* cs,
				const double_array& upDir)
{ // initialize
  assert(0 != _quadrature);
  assert(0 != _db);
  assert(!mesh.isNull());
  assert(0 != cs);

  if (3 != upDir.size())
    throw std::runtime_error("Up direction for surface orientation must be "
			     "a vector with 3 components.");

  // Extract submesh associated with surface
  _boundaryMesh =
    ALE::Selection<ALE::Mesh>::submesh(mesh, mesh->getIntSection(_label));
  if (_boundaryMesh.isNull()) {
    std::ostringstream msg;
    msg << "Could not construct boundary mesh for Neumann traction "
	<< "boundary condition '" << _label << "'.";
    throw std::runtime_error(msg.str());
  } // if

  // check compatibility of quadrature and boundary mesh
  if (_quadrature->cellDim() != _boundaryMesh->getDimension()) {
    std::ostringstream msg;
    msg << "Quadrature is incompatible with cells for Neumann traction "
	<< "boundary condition '" << _label << "'.\n"
	<< "Dimension of boundary mesh: " << _boundaryMesh->getDimension()
	<< ", dimension of quadrature: " << _quadrature->cellDim()
	<< ".";
    throw std::runtime_error(msg.str());
  } // if
  const int numCorners = _quadrature->numBasis();

  // Get 'surface' cells (1 dimension lower than top-level cells)
  const ALE::Obj<Mesh::label_sequence>& cells = _boundaryMesh->heightStratum(1);
  const Mesh::label_sequence::iterator cellsBegin = cells->begin();
  const Mesh::label_sequence::iterator cellsEnd = cells->end();

  // Make sure surface cells are compatible with quadrature.
  const ALE::Obj<sieve_type>& sieve = _boundaryMesh->getSieve();
  assert(!sieve.isNull());
  for (Mesh::label_sequence::iterator c_iter=cellsBegin;
       c_iter != cellsEnd;
       ++c_iter) {
    const int cellNumCorners = sieve->nCone(*c_iter, mesh->depth())->size();
    if (numCorners != cellNumCorners) {
      std::ostringstream msg;
      msg << "Quadrature is incompatible with cell for Neumann traction "
	  << "boundary condition '" << _label << "'.\n"
	  << "Cell " << *c_iter << " has " << cellNumCorners
	  << " vertices but quadrature reference cell has "
	  << numCorners << " vertices.";
      throw std::runtime_error(msg.str());
    } // if
  } // for

  // Create section for traction vector in global coordinates
  const int cellDim = _quadrature->cellDim();
  const int numBasis = _quadrature->numBasis();
  const int numQuadPts = _quadrature->numQuadPts();
  const int spaceDim = cs->spaceDim();
  const int fiberDim = spaceDim * numQuadPts;
  _tractionGlobal = new real_section_type(mesh->comm(), mesh->debug());
  assert(!_tractionGlobal.isNull());
  _tractionGlobal->setFiberDimension(cells, fiberDim);
  mesh->allocate(_tractionGlobal);

  // Containers for orientation information
  const int orientationSize = spaceDim * spaceDim;
  const feassemble::CellGeometry& cellGeometry = _quadrature->refGeometry();
  const int jacobianSize = (cellDim > 0) ? spaceDim * cellDim : 1;
  double_array jacobian(jacobianSize);
  double jacobianDet = 0;
  double_array orientation(orientationSize);
  double_array cellVertices(numBasis*spaceDim);

  // Get mesh coordinates.
  const ALE::Obj<real_section_type>& coordinates =
    mesh->getRealSection("coordinates");

  // open database with traction information
  // NEED TO SET NAMES BASED ON DIMENSION OF BOUNDARY
  // 1-D problem = {'normal-traction'}
  // 2-D problem = {'shear-traction', 'normal-traction'}
  // 3-D problem = {'horiz-shear-traction', 'vert-shear-traction',
  //                'normal-traction'}
  _db->open();
  switch (spaceDim)
    { // switch
    case 1 : {
      const char* valueNames[] = {"normal-traction"};
      _db->queryVals(valueNames, 1);
      break;
    } // case 1
    case 2 : {
      const char* valueNames[] = {"shear-traction", "normal-traction"};
      _db->queryVals(valueNames, 2);
      break;
    } // case 2
    case 3 : {
      const char* valueNames[] = {"horiz-shear-traction", "vert-shear-traction",
				  "normal-traction"};
      _db->queryVals(valueNames, 3);
      break;
    } // case 3
    default :
      assert(0);
    } // switch

  // Containers for database query results and quadrature coordinates in
  // reference geometry.
  double_array tractionDataLocal(spaceDim);
  const double_array& quadPtsRef = _quadrature->quadPtsRef();
  double_array quadPtRef(spaceDim);

  // Container for cell tractions rotated to global coordinates.
  double_array cellTractionsGlobal(fiberDim);

  // Loop over cells in boundary mesh, compute orientations, and then
  // compute corresponding traction vector in global coordinates
  // (store values in _tractionGlobal).
  for(Mesh::label_sequence::iterator c_iter = cellsBegin;
      c_iter != cellsEnd;
      ++c_iter) {
    _quadrature->computeGeometry(_boundaryMesh, coordinates, *c_iter);
    mesh->restrict(coordinates, *c_iter, &cellVertices[0], cellVertices.size());
    const double_array& quadPts = _quadrature->quadPts();

    for(int iQuad = 0, index=0; iQuad < numQuadPts; ++iQuad, index+=spaceDim) {
      // Get traction vector in local coordinate system at quadrature point
      const int err = _db->query(&tractionDataLocal[0], spaceDim,
				 &quadPts[index], spaceDim, cs);
      if (err) {
	std::ostringstream msg;
	msg << "Could not find traction values at \n"
	    << "(";
	for (int i=0; i < spaceDim; ++i)
	  msg << " " << quadPts[index+spaceDim];
	msg << ") for traction boundary condition " << _label << "\n"
	    << "using spatial database " << _db->label() << ".";
	throw std::runtime_error(msg.str());
      } // if

      // Compute Jacobian and determinant at quadrature point, then get
      // orientation.
      memcpy(&quadPtRef[0], &quadPtsRef[index], spaceDim*sizeof(double));
      cellGeometry.jacobian(&jacobian, &jacobianDet, cellVertices, quadPtRef);
      cellGeometry.orientation(&orientation, jacobian, jacobianDet, upDir);

      // Rotate traction vector from local coordinate system to global
      // coordinate system
      cellTractionsGlobal = 0.0;
      for(int iDim = 0; iDim < spaceDim; ++iDim) {
	for(int jDim = 0; jDim < spaceDim; ++jDim)
	  cellTractionsGlobal[iDim] +=
	    orientation[iDim*spaceDim+jDim] * tractionDataLocal[jDim];
      } // for

      // Update tractionGlobal
      mesh->update(_tractionGlobal, *c_iter, &cellTractionsGlobal[0]);
    } // for
  } // for

  _db->close();
} // initialize

// ----------------------------------------------------------------------
// Integrate contributions to residual term (r) for operator.
void
pylith::bc::Neumann::integrateResidual(
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
  const Mesh::label_sequence::iterator cellsBegin = cells->begin();
  const Mesh::label_sequence::iterator cellsEnd = cells->end();

  // Get sections
  const ALE::Obj<real_section_type>& coordinates = 
    mesh->getRealSection("coordinates");
  assert(!coordinates.isNull());

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
  double_array tractionCell(numQuadPts*spaceDim);

  // Loop over faces and integrate contribution from each face
  for (Mesh::label_sequence::iterator c_iter=cellsBegin;
       c_iter != cellsEnd;
       ++c_iter) {
    // Compute geometry information for current cell
    _quadrature->computeGeometry(_boundaryMesh, coordinates, *c_iter);

    // Reset element vector to zero
    _resetCellVector();

    // Restrict tractions to cell
    _boundaryMesh->restrict(_tractionGlobal, *c_iter, 
			    &tractionCell[0], tractionCell.size());

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
            _cellVector[iBasis*spaceDim+iDim] += 
	      tractionCell[iQuad*spaceDim+iDim] * valIJ;
        } // for
      } // for
    } // for
    err = PetscLogFlops(numQuadPts*(1+numBasis*(1+numBasis*(1+2*spaceDim))));
  } // for
} // integrateResidual

// ----------------------------------------------------------------------
// Integrate contributions to Jacobian matrix (A) associated with
void
pylith::bc::Neumann::integrateJacobian(PetscMat* mat,
				       const double t,
				       topology::FieldsManager* const fields,
				       const ALE::Obj<Mesh>& mesh)
{ // integrateJacobian
  throw std::logic_error("Not implemented.");
} // integrateJacobian

// ----------------------------------------------------------------------
// Verify configuration is acceptable.
void
pylith::bc::Neumann::verifyConfiguration(const ALE::Obj<Mesh>& mesh)
{ // verifyConfiguration
  throw std::logic_error("Not implemented.");
} // verifyConfiguration


// End of file 
