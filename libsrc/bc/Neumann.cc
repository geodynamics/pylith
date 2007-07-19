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

#include "spatialdata/spatialdb/SpatialDB.hh" // USES SpatialDB
#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys

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
				const double_array& upDir);
{ // initialize
  assert(0 != _quadrature);
  assert(0 != _db);
  assert(!mesh.isNull());
  assert(0 != cs);

  if (3 != upDir.size())
    throw std::runtime_error("Up direction for surface orientation must be "
			     "a vector with 3 components.");

  // Extract submesh associated with surface
  const ALE::Obj<ALE::Mesh> boundaryMesh =
    ALE::Selection<ALE::Mesh>::submesh(mesh, mesh->getIntSection(_label));
  if (boundaryMesh.isNull()) {
    std::ostringstream msg;
    msg << "Could not construct boundary mesh for Neumann traction "
	<< "boundary condition '" << _label << "'.";
    throw std::runtime_error(msg.str());
  } // if

  // Get 'surface' cells (1 dimension lower than top-level cells)
  const ALE::Obj<sieve_type>& sieve = boundaryMesh->getSieve();
  const ALE::Obj<Mesh::label_sequence>& cells = boundaryMesh->heightStratum(1);
  const int numCells = cells->size();

  // Create section for traction vector in global coordinates
  const int spaceDim = cs->spaceDim();
  const int cellDim = _quadrature->cellDim();
  const int numQuadPts = _quadrature->numQuadPts();
  const int fiberDim = spaceDim * numQuadPts;
  _tractionGlobal = new real_section_type(mesh->comm(), mesh->debug());
  assert(!_tractionGlobal.isNull());
  _tractionGlobal->setFiberDimension(cells, fiberDim);
  _tractionGlobal->allocate();

  const int numBasis = _quadrature->numBasis();
  const feassemble::CellGeometry& cellGeometry = _quadrature->refGeometry();
  const double_array& verticesRef = _quadrature->vertices();
  const int jacobianSize = (cellDim > 0) ? spaceDim * cellDim : 1;
  double_array jacobian(jacobianSize);
  double jacobianDet = 0;
  double_array cellVertices(numBasis*spaceDim);

  // Set up orientation information
  const int orientationSize = spaceDim * spaceDim;
  double_array orientation(orientationSize);
  const ALE::Obj<real_section_type>& coordinates =
    mesh->getRealSection("coordinates");

  // setup database with traction information
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

  // Loop over cells in boundary mesh, compute orientations, and then
  // compute corresponding traction vector in global coordinates
  // (store values in _tractionGlobal).
  const Mesh::label_sequence::iterator cellsBegin = cells->begin();
  const Mesh::label_sequence::iterator cellsEnd = cells->end();
  double_array tractionDataLocal(spaceDim);
  const feassemble::CellGeometry& cellGeometry = _quadrature->refGeometry();
  const double_array& quadPtsRef = _quadrature->quadPtsRef();
  double_array rotatedTraction[spaceDim*numQuadPts];

  for(Mesh::label_sequence::iterator c_iter = cells->begin();
      c_iter != cells->end();
      ++c_iter) {
    _quadrature->computeGeometry(mesh, coordinates, *c_iter);
    mesh->restrict(coordinates, *c_iter, &cellVertices[0], cellVertices.size());
    const double_array& quadPts = _quadrature->quadPts();
    for(int iQuad = 0, iQ=iQuad*spaceDim; iQuad < numQuadPts; ++iQuad) {
      // Get traction vector in local coordinate system at quadrature point
      int err = _db->query(&tractionDataLocal[0], spaceDim, quadPts[iQ + iQuad],
			   spaceDim, cs);

      // Compute Jacobian and determinant at quadrature point
      double_array quadPt(&quadPtsRef[iQuad*cellDim], cellDim);
      cellGeometry.jacobian(&jacobian, &jacobianDet, cellVertices, quadPt);
      
      // Get orientation of boundary at quadrature point
      cellGeometry.orientation(&orientation, jacobian, jacobianDet, upDir);

      // Normalize orientation
      double sum;
      double_array quadPtOrient(orientationSize);
      for(int iDim = 0; iDim < spaceDim; ++iDim) {
	sum = 0.0;
	for(int jDim = 0, index=iDim*spaceDim; jDim < spaceDim; ++jDim)
	  sum += pow(orientation[index+jDim], 2);
	sum =sqrt(sum);
	assert(sum > 0.0);
	for (int jDim=0, index=iDim*spaceDim; jDim < spaceDim; ++jDim)
	  quadPtOrient[index+jDim] = orientation[index+jDim] / sum;
      } // for

      // Rotate traction vector from local coordinate system to global
      // coordinate system
      for(int iDim = 0; iDim < spaceDim; ++iDim) {
	sum = 0.0;
	for(int jDim = 0, index=iDim*spaceDim; jDim < spaceDim; ++jDim)
	  sum += orientation[index+jDim] * tractionDataLocal[jDim];
	rotatedTraction[iDim+iQ] = sum;
      } // for

      // Update tractionGlobal
      _tractionGlobal->updateAddPoint(*c_iter, &rotatedTraction[0]);
    } // for
  } // for

  _db->close();
} // initialize

// ----------------------------------------------------------------------
// Integrate contributions to residual term (r) for operator.
void
pylith::bc::Neumann::integrateResidual(
				  const ALE::Obj<real_section_type>& residual,
				  const ALE::Obj<Mesh>& mesh)
{ // integrateResidual
  throw std::logic_error("Not implemented.");
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
