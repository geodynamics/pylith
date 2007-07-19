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

#include "spatialdata/spatialdb/SpatialDB.hh" // USES SpatialDB
#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys

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
				const double_array& upDir);
{ // initialize
  assert(0 != _quadrature);
  assert(0 != _db);
  assert(!mesh.isNull());
  assert(0 != cs);

  if (3 != upDir.size())
    throw std::runtime_error("Up direction for boundary orientation must be "
			     "a vector with 3 components.");

  // Extract submesh associated with boundary
  const ALE::Obj<ALE::Mesh> boundaryMesh =
    ALE::Selection<ALE::Mesh>::submesh(mesh, mesh->getIntSection(_label));
  if (boundaryMesh.isNull()) {
    std::ostringstream msg;
    msg << "Could not construct boundary mesh for absorbing boundary "
	<< "condition '" << _label << "'.";
    throw std::runtime_error(msg.str());
  } // if

  // Get 'surface' cells (1 dimension lower than top-level cells)
  const ALE::Obj<sieve_type>& sieve = boundaryMesh->getSieve();
  const ALE::Obj<Mesh::label_sequence>& cells = boundaryMesh->heightStratum(1);

  // Get damping constants at each quadrature point
  const int spaceDim = cs->spaceDim();
  const int numQuadPts = _quadrature->numQuadPts();
  const int fiberDim = spaceDim * numQuadPts;
  _dampingConsts = new real_section_type(mesh->comm(), mesh->debug());
  assert(!_dampingConsts.isNull());
  _dampingConsts->setFiberDimension(cells, fiberDim);
  _dampingConsts->allocate();

  const int numBasis = _quadrature->numBasis();
  const feassemble::CellGeometry& cellGeometry = _quadrature->refGeometry();
  const double_array& verticesRef = _quadrature->vertices();
  const int jacobianSize = (cellDim > 0) ? spaceDim * cellDim : 1;
  double_array jacobian(jacobianSize);
  double jacobianDet = 0;
  double_array cellVertices(numBasis*spaceDim);

  // Set up orientation information
  const int orientationSize = spaceDim * spaceDim;
  double orientation[orientationSize];
  const ALE::Obj<real_section_type>& coordinates =
    mesh->getRealSection("coordinates");

  // open database with material property information
  _db->open();
  const char* valuesNames = 
    (cellDim > 0) ? { "vp", "vs", "density" } : { "vp", "density" };
  const int numValues = (cellDim > 0) ? 3 : 2;
  _db->queryVals((const char**) valueNames, numValues);

  const Mesh::label_sequence::iterator cellsBegin = cells->begin();
  const Mesh::label_sequence::iterator cellsEnd = cells->end();
  for(Mesh::label_sequence::iterator c_iter = cells->begin();
      c_iter != cells->end();
      ++c_iter) {
    _quadrature->computeGeometry(mesh, coordinates, *c_iter);
    for(int iQuad = 0; iQuad < numQuadPts; ++iQuad) {
    } // for
  } // for

  _db->close();
} // initialize

// ----------------------------------------------------------------------
// Integrate contributions to residual term (r) for operator.
void
pylith::bc::AbsorbingDampers::integrateResidual(
				  const ALE::Obj<real_section_type>& residual,
				  const ALE::Obj<Mesh>& mesh)
{ // integrateResidual
  throw std::logic_error("Not implemented.");
} // integrateResidual

// ----------------------------------------------------------------------
// Integrate contributions to Jacobian matrix (A) associated with
void
pylith::bc::AbsorbingDampers::integrateJacobian(
					 PetscMat* mat,
					 const double t,
					 topology::FieldsManager* const fields,
					 const ALE::Obj<Mesh>& mesh)
{ // integrateJacobian
  throw std::logic_error("Not implemented.");
} // integrateJacobian

// ----------------------------------------------------------------------
// Verify configuration is acceptable.
void
pylith::bc::AbsorbingDampers::verifyConfiguration(const ALE::Obj<Mesh>& mesh)
{ // verifyConfiguration
  throw std::logic_error("Not implemented.");
} // verifyConfiguration


// End of file 
